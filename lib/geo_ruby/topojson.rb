require 'json'
require 'geo_ruby/simple_features'
require 'geo_ruby/geojson'

module GeoRuby
  class TopoJSONError < StandardError
  end

  class TopoJSONExporter
    DEFAULT_EXTENT = 10_000

    GeometryRecord = Struct.new(:kind, :coordinates, :geometries, :properties, :id, keyword_init: true)

    PathRecord = Struct.new(:points, :closed, keyword_init: true)

    def self.export(input, object_name: 'data', extent: DEFAULT_EXTENT, transform: nil)
      new(extent: extent, transform: transform).export(input, object_name: object_name)
    end

    def initialize(extent: DEFAULT_EXTENT, transform: nil)
      raise TopoJSONError, 'extent must be at least 2' unless extent.is_a?(Integer) && extent >= 2

      @extent = extent
      @transform = validate_transform(transform)
      @diagnostics = []
      @arc_index = {}
      @absolute_arcs = []
      @collection_input = false
    end

    def export(input, object_name: 'data')
      features = normalize_input(input)
      transform = @transform || build_transform(features)
      @transform = transform
      records = features.map.with_index do |feature, index|
        record = prepare_feature(feature)
        properties = record.properties
        id = record.id
        quantized = quantize_record(record, [index])
        quantized.properties = properties if quantized
        quantized.id = id if quantized
        quantized
      end

      paths = []
      records.each { |record| collect_paths(record, paths) }
      split_points = collect_split_points(paths)
      paths.each { |path| split_path!(path, split_points) }

      object = if @collection_input
        { 'type' => 'GeometryCollection', 'geometries' => records.map { |record| convert_record(record) } }
      else
        convert_record(records.first)
      end

      topology = stringify_keys({
        'type' => 'Topology',
        'transform' => serialize_transform(transform),
        'objects' => { object_name => object },
        'arcs' => encoded_arcs
      })
      topology['diagnostics'] = stringify_keys(@diagnostics) unless @diagnostics.empty?
      topology
    end

    private

    def normalize_input(input)
      return [] if input.nil?

      if input.is_a?(GeoJSONFeatureCollection)
        @collection_input = true
        input.features.to_a
      elsif input.is_a?(Array)
        @collection_input = true
        input.map do |item|
          item.is_a?(Array) ? GeoJSONFeature.new(item[0], item[1] || {}, item[2]) : item
        end
      else
        @collection_input = false
        [input.is_a?(GeoJSONFeature) ? input : GeoJSONFeature.new(input, {})]
      end
    end

    def prepare_feature(feature)
      if feature.is_a?(GeoJSONFeature)
        record = prepare_geometry(feature.geometry, [])
        record.properties = feature.properties || {}
        record.id = feature.id unless feature.id.nil?
        record
      else
        prepare_geometry(feature, [])
      end
    end

    def prepare_geometry(geometry, path)
      return null_record if geometry.nil?

      case geometry
      when SimpleFeatures::Point
        coordinates = prepare_coordinates(geometry.to_coordinates, path + ['coordinates'])
        coordinates.nil? ? null_record : GeometryRecord.new(kind: 'Point', coordinates: coordinates[0, 2])
      when SimpleFeatures::MultiPoint
        coordinates = geometry.geometries.map.with_index do |point, index|
          prepared = prepare_coordinates(point.to_coordinates, path + [index])
          break nil if prepared.nil?
          prepared[0, 2]
        end
        return null_record if coordinates.nil?
        GeometryRecord.new(kind: 'MultiPoint', coordinates: coordinates)
      when SimpleFeatures::LineString
        prepare_line(geometry.points, false)
      when SimpleFeatures::MultiLineString
        lines = geometry.geometries.map do |line|
          record = prepare_line(line.points, false)
          break nil if null_record?(record)
          record.coordinates
        end
        return null_record if lines.nil?
        GeometryRecord.new(kind: 'MultiLineString', coordinates: lines)
      when SimpleFeatures::Polygon
        prepare_polygon(geometry.rings, path)
      when SimpleFeatures::MultiPolygon
        polygons = geometry.geometries.map.with_index do |polygon, index|
          record = prepare_polygon(polygon.rings, path + [index])
          break nil if null_record?(record)
          record.coordinates
        end
        return null_record if polygons.nil?
        GeometryRecord.new(kind: 'MultiPolygon', coordinates: polygons)
      when SimpleFeatures::GeometryCollection
        members = geometry.geometries.map.with_index { |member, index| prepare_geometry(member, path + [index]) }
        GeometryRecord.new(kind: 'GeometryCollection', geometries: members)
      else
        diagnostic(path, 'unsupported_geometry', "#{geometry.class} cannot be represented as TopoJSON")
        null_record
      end
    end

    def prepare_line(points, closed)
      return GeometryRecord.new(kind: 'LineString', coordinates: []) if !closed && points.empty?

      coordinates = points.map { |point| point.to_coordinates[0, 2] }
      GeometryRecord.new(kind: closed ? 'Polygon' : 'LineString', coordinates: [PathRecord.new(points: coordinates, closed: closed)])
    end

    def prepare_polygon(rings, path)
      if rings.empty?
        diagnostic(path, 'empty_polygon', 'polygon has no rings')
        return null_record
      end

      paths = rings.map do |ring|
        record = prepare_line(ring.points, true)
        break nil if null_record?(record)
        record.coordinates.first
      end
      return null_record if paths.nil?
      GeometryRecord.new(kind: 'Polygon', coordinates: paths)
    end

    def prepare_coordinates(coordinates, path)
      if !coordinates.is_a?(Array) || coordinates.length < 2 || !finite_number?(coordinates[0]) || !finite_number?(coordinates[1])
        diagnostic(path, 'invalid_coordinates', 'x and y must be finite numbers')
        return nil
      end
      coordinates
    end

    def finite_number?(value)
      value.is_a?(Numeric) && value.finite?
    end

    def null_record
      GeometryRecord.new(kind: nil)
    end

    def null_record?(record)
      record.kind.nil?
    end

    def build_transform(features)
      min_x = min_y = Float::INFINITY
      max_x = max_y = -Float::INFINITY

      each_geometry(features) do |geometry|
        each_point(geometry) do |point|
          next unless point.x.finite? && point.y.finite?

          x = point.x.to_f
          y = point.y.to_f
          min_x = [min_x, x].min
          max_x = [max_x, x].max
          min_y = [min_y, y].min
          max_y = [max_y, y].max
        end
      end

      if min_x.infinite?
        return { scale: [1.0, 1.0], translate: [0.0, 0.0] }
      end

      width = [max_x - min_x, 0.0].max
      height = [max_y - min_y, 0.0].max
      {
        scale: [
          width.zero? ? 1.0 : width / (@extent - 1),
          height.zero? ? 1.0 : height / (@extent - 1)
        ],
        translate: [min_x, min_y]
      }
    end

    def each_geometry(features)
      features.each do |feature|
        yield(feature.is_a?(GeoJSONFeature) ? feature.geometry : feature)
      end
    end

    def each_point(geometry)
      return if geometry.nil?

      case geometry
      when SimpleFeatures::Point
        yield geometry
      when SimpleFeatures::LineString
        geometry.points.each { |point| yield point }
      when SimpleFeatures::Polygon
        geometry.rings.each { |ring| ring.points.each { |point| yield point } }
      when SimpleFeatures::GeometryCollection
        geometry.geometries.each { |member| each_point(member) { |point| yield point } }
      end
    end

    def validate_transform(transform)
      return nil if transform.nil?

      unless transform.is_a?(Hash) && transform[:scale].is_a?(Array) && transform[:translate].is_a?(Array) &&
             transform[:scale].length == 2 && transform[:translate].length == 2 &&
             transform[:scale].all? { |value| finite_number?(value) && value > 0 } &&
             transform[:translate].all? { |value| finite_number?(value) }
        raise TopoJSONError, 'transform must contain positive scale and finite translate values'
      end
      transform
    end

    def quantize_record(record, path)
      return record if null_record?(record)

      case record.kind
      when 'Point'
        record.coordinates = quantize_coordinates(record.coordinates, path + ['coordinates'])
        return nullify(record, path, 'invalid_coordinates', 'point has unrepresentable coordinates') if record.coordinates.nil?
      when 'MultiPoint'
        record.coordinates = record.coordinates.map.with_index do |coordinates, index|
          quantized = quantize_coordinates(coordinates, path + [index])
          break nil if quantized.nil?
          quantized
        end
        return nullify(record, path, 'invalid_coordinates', 'multi point contains unrepresentable coordinates') if record.coordinates.nil?
      when 'LineString'
        return record if record.coordinates == []
        record.coordinates = record.coordinates.map do |path_record|
          quantize_path_record(path_record, path)
        end
        non_empty = record.coordinates.compact
        return nullify(record, path, 'degenerate_line_string', 'line string has no non-degenerate edges') if non_empty.empty?
        record.coordinates = non_empty
      when 'Polygon'
        record.coordinates = record.coordinates.map.with_index do |path_record, index|
          quantize_path_record(path_record, path + [index])
        end.compact
        return nullify(record, path, 'degenerate_polygon', 'polygon has no non-degenerate rings') if record.coordinates.empty?
      when 'MultiLineString'
        record.coordinates = record.coordinates.map.with_index do |line, line_index|
          line.map { |path_record| quantize_path_record(path_record, path + [line_index]) }.compact
        end
        return nullify(record, path, 'degenerate_line_string', 'multi line string has no non-degenerate lines') if record.coordinates.any?(&:empty?)
      when 'MultiPolygon'
        record.coordinates = record.coordinates.map.with_index do |polygon, polygon_index|
          rings = polygon.map.with_index do |path_record, ring_index|
            quantize_path_record(path_record, path + [polygon_index, ring_index])
          end.compact
          break nil if rings.empty?
          rings
        end
        return nullify(record, path, 'degenerate_polygon', 'multi polygon contains an empty polygon') if record.coordinates.nil?
      when 'GeometryCollection'
        record.geometries = record.geometries.map.with_index do |member, index|
          quantize_record(member, path + [index])
        end
      end
      record
    end

    def quantize_coordinates(coordinates, path)
      unless coordinates.is_a?(Array) && coordinates.length >= 2 &&
           finite_number?(coordinates[0]) && finite_number?(coordinates[1])
        diagnostic(path, 'invalid_coordinates', 'x and y must be finite numbers')
        return nil
      end

      x = quantize_value(coordinates[0], @transform[:scale][0], @transform[:translate][0])
      y = quantize_value(coordinates[1], @transform[:scale][1], @transform[:translate][1])
      if x < 0 || x > @extent - 1 || y < 0 || y > @extent - 1
        diagnostic(path, 'coordinate_outside_extent', "coordinate is outside quantized extent 0..#{@extent - 1}")
        return nil
      end
      [x, y]
    end

    def quantize_value(value, scale, translate)
      if scale.zero?
        0
      else
        ((value.to_f - translate) / scale).round
      end
    end

    def quantize_path_record(path_record, path)
      return PathRecord.new(points: [], closed: false) if !path_record.closed && path_record.points.empty?

      points = path_record.points.map do |coordinates|
        quantize_coordinates(coordinates, path + ['coordinates'])
      end
      return if points.include?(nil)

      if path_record.closed && points.first != points.last
        diagnostic(path, 'unclosed_ring', 'ring is not closed')
        return nil
      end
      points = points[0...-1] if path_record.closed
      points = [points.first] + points.each_cons(2).map do |previous, current|
        current unless current == previous
      end
      points.compact!

      minimum_points = path_record.closed ? 3 : 2
      if points.uniq.length < minimum_points
        diagnostic(path, 'degenerate_geometry', 'quantization removed all distinct edges')
        return nil
      end

      PathRecord.new(points: points, closed: path_record.closed)
    end

    def nullify(record, path, code, message)
      diagnostic(path, code, message)
      properties = record.properties
      id = record.id
      record.kind = nil
      record.coordinates = nil
      record.geometries = nil
      record.properties = properties
      record.id = id
      record
    end

    def collect_paths(record, paths)
      return if null_record?(record)

      case record.kind
      when 'LineString'
        paths.concat(record.coordinates)
      when 'Polygon'
        paths.concat(record.coordinates)
      when 'MultiLineString'
        record.coordinates.each { |line| paths.concat(line) }
      when 'MultiPolygon'
        record.coordinates.each { |polygon| paths.concat(polygon) }
      when 'GeometryCollection'
        record.geometries.each { |member| collect_paths(member, paths) }
      end
    end

    def collect_split_points(paths)
      all_points = Set.new
      paths.each do |path|
        path.points.each do |point|
          point = point.dup
          all_points.add(point) unless all_points.include?(point)
        end
      end
      all_points
    end

    def cross_product(start, endpoint, point)
      (endpoint.first - start.first) * (point.last - start.last) -
        (endpoint.last - start.last) * (point.first - start.first)
    end

    def split_path!(path, split_points)
      split = []
      path.points.each_cons(2) do |start, endpoint|
        next if compare_points(start, endpoint).zero?
        start = [start.first, start.last]
        endpoint = [endpoint.first, endpoint.last]
        split << start if split.empty? || compare_points(split.last, start) != 0
        between = split_points.select do |point|
          strictly_inside_segment?(start, endpoint, point)
        end
        direction = start.first != endpoint.first ? start.first <=> endpoint.first : start.last <=> endpoint.last
        between.sort_by { |point| start.first != endpoint.first ? point.first * direction : point.last * direction }.each do |point|
          split << [point.first, point.last]
        end
        split << endpoint
      end
      path.points = split
    end

    def point_on_segment?(start, endpoint, point)
      x_range = [start.first, endpoint.first].minmax
      y_range = [start.last, endpoint.last].minmax
      point.first.between?(*x_range) && point.last.between?(*y_range) &&
        cross_product(start, endpoint, point).zero?
    end

    def strictly_inside_segment?(start, endpoint, point)
      point != start && point != endpoint && point_on_segment?(start, endpoint, point)
    end

    def convert_record(record)
      return null_geometry(record) if null_record?(record)

      case record.kind
      when 'Point'
        base_record(record, 'coordinates', record.coordinates)
      when 'MultiPoint'
        base_record(record, 'coordinates', record.coordinates)
      when 'LineString'
        arcs = record.coordinates.map { |path| path_arcs(path) }
        base_record(record, 'arcs', arcs)
      when 'Polygon'
        arcs = record.coordinates.map { |path| path_arcs(path) }
        base_record(record, 'arcs', arcs)
      when 'MultiLineString'
        arcs = record.coordinates.map { |line| line.map { |path| path_arcs(path) } }
        base_record(record, 'arcs', arcs)
      when 'MultiPolygon'
        arcs = record.coordinates.map do |polygon|
          polygon.map { |path| path_arcs(path) }
        end
        base_record(record, 'arcs', arcs)
      when 'GeometryCollection'
        {
          'type' => 'GeometryCollection',
          'geometries' => record.geometries.map { |member| convert_record(member) },
          'properties' => record.properties
        }.tap { |object| object['id'] = record.id unless record.id.nil? }
      end
    end

    def base_record(record, coordinate_key, value)
      {
        'type' => record.kind,
        coordinate_key => value,
        'properties' => record.properties
      }.tap { |object| object['id'] = record.id unless record.id.nil? }
    end

    def null_geometry(record)
      {
        'type' => nil,
        'coordinates' => nil,
        'properties' => record.properties
      }.tap { |object| object['id'] = record.id unless record.id.nil? }
    end

    def path_arcs(path)
      ring_points = path.closed ? path.points + [path.points.first] : path.points
      ring_points.each_cons(2).map do |start, endpoint|
        arc_reference(start, endpoint)
      end
    end

    def arc_reference(start, endpoint)
      if compare_points(start, endpoint).negative?
        canonical_start, canonical_end = start, endpoint
      else
        canonical_start, canonical_end = endpoint, start
      end
      reversed = compare_points(start, canonical_start) != 0
      canonical_start = [canonical_start.first, canonical_start.last]
      canonical_end = [canonical_end.first, canonical_end.last]
      key = arc_key(canonical_start, canonical_end)
      index = @arc_index[key]
      if index.nil?
        index = @absolute_arcs.length
        @absolute_arcs << [canonical_start.dup, canonical_end.dup]
        @arc_index[key] = index
      end
      reversed ? ~index : index
    end

    def compare_points(first, second)
      (first.first <=> second.first).nonzero? || first.last <=> second.last
    end

    def arc_key(start, endpoint)
      [start.first, start.last, endpoint.first, endpoint.last]
    end

    def encoded_arcs
      @absolute_arcs.map do |points|
        previous = [0, 0]
        points.map do |point|
          delta = [point.first - previous.first, point.last - previous.last]
          previous = point.dup
          delta
        end
      end
    end

    def serialize_transform(transform)
      {
        'scale' => transform[:scale],
        'translate' => transform[:translate]
      }
    end

    def diagnostic(path, code, message)
      @diagnostics << { 'path' => path, 'code' => code, 'message' => message }
    end

    def stringify_keys(value)
      case value
      when Hash
        value.each_with_object({}) do |(key, item), result|
          result[key.to_s] = stringify_keys(item)
        end
      when Array
        value.map { |item| stringify_keys(item) }
      else
        value
      end
    end
  end
end
