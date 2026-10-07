# frozen_string_literal: true

require 'geo_ruby/simple_features'
require 'geo_ruby/geojson'

begin
  require 'json'
rescue LoadError
  warn "Please install a JSON provider: gem install json"
  raise
end

module GeoRuby
  class TopologyResult
    attr_reader :topology, :diagnostics

    def initialize(topology, diagnostics)
      @topology = topology
      @diagnostics = diagnostics
    end

    def to_json(*arguments)
      @topology.to_json(*arguments)
    end

    def as_json(*arguments)
      @topology.as_json(*arguments)
    end

    alias_method :to_topojson, :to_json
  end

  class TopoJSONExporter
    # Exports GeoRuby geometries and GeoJSON features as a quantized
    # TopoJSON topology. The processing order is deliberate:
    # 1. quantize all coordinates to transform-space integers;
    # 2. remove consecutive points merged by quantization;
    # 3. node and deduplicate directed edges;
    # 4. stitch maximal shared arcs.
    GeometryClass = SimpleFeatures::Geometry
    PointClass = SimpleFeatures::Point
    LineStringClass = SimpleFeatures::LineString
    LinearRingClass = SimpleFeatures::LinearRing
    PolygonClass = SimpleFeatures::Polygon
    MultiPointClass = SimpleFeatures::MultiPoint
    MultiLineStringClass = SimpleFeatures::MultiLineString
    MultiPolygonClass = SimpleFeatures::MultiPolygon
    GeometryCollectionClass = SimpleFeatures::GeometryCollection
    CircleClass = SimpleFeatures::Circle

    def self.export(input, object_name: 'collection', quantization: 10_000)
      new(input, object_name: object_name, quantization: quantization).export
    end

    def initialize(input, object_name: 'collection', quantization: 10_000)
      unless quantization.is_a?(Integer) && quantization >= 2
        raise ArgumentError, 'quantization must be an integer greater than 1'
      end
      @input = input
      @object_name = object_name
      @quantization = quantization
      @diagnostics = []
    end

    def export
      objects = normalize(@input)

      bounds = collect_bounds(objects.values)
      @transform = build_transform(bounds)
      quantize_objects(objects)
      topology = build_topology(objects, @transform)
      TopologyResult.new(topology, @diagnostics)
    end

    private

    def normalize(input)
      if input.is_a?(GeoJSONFeatureCollection)
        { @object_name => {
          'type' => 'GeometryCollection',
          'geometries' => input.features.map { |feature| prepare_feature(feature) }
        } }
      elsif input.is_a?(GeoJSONFeature)
        { @object_name => prepare_feature(input) }
      elsif input.is_a?(Array)
        { @object_name => {
          'type' => 'GeometryCollection',
          'geometries' => input.map { |item| wrap_feature(item) }
        } }
      elsif input.is_a?(Hash)
        input.each_with_object({}) do |(name, item), objects|
          objects[name.to_s] = wrap_feature(item)
        end
      elsif input.is_a?(GeometryClass)
        { @object_name => prepare_geometry(input) }
      elsif input.nil?
        { @object_name => null_object }
      else
        raise ArgumentError, "unsupported TopoJSON input: #{input.class}"
      end
    end

    def wrap_feature(item)
      if item.is_a?(GeoJSONFeature)
        prepare_feature(item)
      elsif item.nil?
        null_object
      else
        prepare_geometry(item)
      end
    end

    def prepare_feature(feature)
      object = if feature.geometry.nil?
                 null_object
               else
                 prepare_geometry(feature.geometry)
               end
      object['id'] = feature.id unless feature.id.nil?
      object['properties'] = feature.properties
      object
    end

    def prepare_geometry(geometry)
      if geometry.is_a?(PointClass)
        { 'type' => 'Point', 'coordinates' => coordinates_for(geometry) }
      elsif geometry.is_a?(MultiPointClass)
        { 'type' => 'MultiPoint', 'coordinates' => coordinates_for(geometry) }
      elsif geometry.is_a?(LinearRingClass) || geometry.is_a?(LineStringClass)
        { 'type' => 'LineString', 'coordinates' => coordinates_for(geometry) }
      elsif geometry.is_a?(PolygonClass)
        { 'type' => 'Polygon',
          'arcs' => geometry.rings.map { |ring| coordinates_for(ring) } }
      elsif geometry.is_a?(MultiLineStringClass)
        lines = geometry.geometries.map { |line| coordinates_for(line) }
        { 'type' => 'MultiLineString', 'arcs' => lines }
      elsif geometry.is_a?(MultiPolygonClass)
        polygons = geometry.geometries.map do |polygon|
          polygon.rings.map { |ring| coordinates_for(ring) }
        end
        { 'type' => 'MultiPolygon', 'arcs' => polygons }
      elsif geometry.is_a?(GeometryCollectionClass)
        children = geometry.geometries.map { |child| prepare_geometry(child) }
        { 'type' => 'GeometryCollection', 'geometries' => children }
      elsif geometry.is_a?(CircleClass)
        object = null_object
        object['_diagnostic'] = { reason: :unsupported_geometry, geometry: 'Circle' }
        object
      elsif geometry.nil?
        null_object
      else
        object = null_object
        object['_diagnostic'] = {
          reason: :unsupported_geometry,
          geometry: geometry.class.name
        }
        object
      end
    end

    def coordinates_for(geometry)
      geometry.to_coordinates
    end

    def null_object
      { 'type' => nil }
    end

    def finite_coordinates?(*coordinates)
      coordinates.all? { |coordinate| coordinate.is_a?(Numeric) && coordinate.finite? }
    end

    def nullify(object)
      object['type'] = nil
      object.delete('coordinates')
      object.delete('arcs')
      object.delete('geometries')
      object
    end

    def collect_bounds(objects)
      points = []
      objects.each { |object| collect_positions(object, points) }
      return nil if points.empty?
      xs = points.map(&:first)
      ys = points.map(&:last)
      [xs.min, ys.min, xs.max, ys.max]
    end

    def collect_positions(object, points)
      case object['type']
      when 'Point' then points << object['coordinates']
      when 'MultiPoint' then points.concat(object['coordinates'])
      when 'LineString' then points.concat(object['coordinates'])
      when 'MultiLineString' then object['arcs'].each { |line| points.concat(line) }
      when 'Polygon' then object['arcs'].each { |ring| points.concat(ring) }
      when 'MultiPolygon'
        object['arcs'].each do |polygon|
          polygon.each { |ring| points.concat(ring) }
        end
      when 'GeometryCollection'
        object['geometries'].each { |child| collect_positions(child, points) }
      end
    end

    def build_transform(bounds)
      if bounds.nil?
        { 'scale' => [1, 1], 'translate' => [0, 0] }
      else
        min_x, min_y, max_x, max_y = bounds
        width = max_x - min_x
        height = max_y - min_y
        scale_x = width.zero? ? 1 : width.to_f / (@quantization - 1)
        scale_y = height.zero? ? 1 : height.to_f / (@quantization - 1)
        { 'scale' => [scale_x, scale_y], 'translate' => [min_x, min_y] }
      end
    end

    def quantize_objects(objects)
      objects.each do |name, object|
        quantize_object(object, base_context(object, name: name))
      end
    end

    def quantize_object(object, context = nil)
      context ||= base_context(object)
      if object['_diagnostic']
        diagnostic(object['_diagnostic'][:reason], context.merge(object['_diagnostic']))
        return
      end

      case object['type']
      when 'Point' then quantize_point(object, context)
      when 'MultiPoint' then quantize_multi_point(object, context)
      when 'LineString' then quantize_line(object, context)
      when 'MultiLineString' then quantize_multi_line(object, context)
      when 'Polygon' then quantize_polygon(object, context)
      when 'MultiPolygon' then quantize_multi_polygon(object, context)
      when 'GeometryCollection'
        object['geometries'].each_with_index do |child, index|
          quantize_object(child, child_context(child, context, index))
        end
      end
    end

    def quantize_point(object, context)
      valid = object['coordinates'].all? do |coordinate|
        coordinate.is_a?(Numeric) && coordinate.finite?
      end
      unless valid
        diagnostic(:invalid_coordinates, context.merge(geometry: 'Point'))
        nullify(object)
        return
      end
      object['coordinates'] = quantize_position(object['coordinates'])
    end

    def quantize_multi_point(object, context)
      unless object['coordinates'].all? { |position| finite_coordinates?(*position) }
        diagnostic(:invalid_coordinates, context.merge(geometry: 'MultiPoint'))
        nullify(object)
        return
      end
      object['coordinates'].map! { |position| quantize_position(position) }
    end

    def quantize_line(object, context, geometry: 'LineString')
      path = object.delete('coordinates')
      positions, error = quantize_path(path)
      if error
        diagnostic(error || :degenerate_line, context.merge(geometry: geometry))
        nullify(object)
      elsif positions.one?
        diagnostic(:degenerate_line, context.merge(geometry: geometry))
        nullify(object)
      else
        object['arcs'] = positions
      end
    end

    def quantize_multi_line(object, context)
      lines = object['arcs'].map do |line|
        positions, error = quantize_path(line)
        if error
          diagnostic(error, context.merge(geometry: 'MultiLineString'))
          nil
        elsif positions.length < 2
          diagnostic(:degenerate_line, context.merge(geometry: geometry))
          nil
        else
          positions
        end
      end.compact
      object['arcs'] = lines
    end

    def quantize_polygon(object, context)
      if object['arcs'].empty?
        object['arcs'] = []
        return
      end
      outer_valid = false
      rings = object['arcs'].map.with_index do |ring, ring_index|
        quantized_ring, error = quantize_ring(ring)
        if error
          diagnostic(error, context.merge(geometry: 'Polygon', ring: ring_index))
          nil
        elsif ring_invalid?(quantized_ring, ring_index, context)
          nil
        else
          outer_valid = true if ring_index.zero?
          quantized_ring
        end
      end.compact
      if outer_valid
        object['arcs'] = rings
      else
        nullify(object)
      end
    end

    def quantize_multi_polygon(object, context)
      polygons = object['arcs'].map.with_index do |rings, polygon_index|
        outer_valid = false
        quantized_rings = rings.map.with_index do |ring, ring_index|
          quantized_ring, error = quantize_ring(ring)
          if error
            diagnostic(
              error,
              context.merge(
                geometry: 'MultiPolygon',
                polygon: polygon_index,
                ring: ring_index
              )
            )
            nil
          elsif ring_invalid?(quantized_ring, ring_index, context.merge(
            geometry: 'Polygon', polygon: polygon_index
          ))
            nil
          else
            outer_valid = true if ring_index.zero?
            quantized_ring
          end
        end.compact
        outer_valid || quantized_rings.empty? ? quantized_rings : nil
      end
      object['arcs'] = polygons.compact
    end

    def quantize_ring(ring)
      positions, error = quantize_path(ring)
      return [[], error] if error
      return [positions, nil] if positions.empty?
      positions << positions.first unless positions.last == positions.first
      positions = positions.chunk { |position| position }.map(&:first)
      if positions.length > 1 && positions.last != positions.first
        positions << positions.first
      end
      if positions.length < 4 || ring_area(positions).zero?
        return [[], :degenerate_ring]
      end
      [positions, nil]
    end

    def quantize_path(path)
      unless path.all? { |position| finite_coordinates?(*position) }
        return [[], :invalid_coordinates]
      end
      positions = path.map { |position| quantize_position(position) }
      [positions.chunk { |position| position }.map(&:first), nil]
    end

    def quantize_position(position)
      scale = @transform['scale']
      translate = @transform['translate']
      quantized = [
        quantize_coordinate(position[0], translate[0], scale[0]),
        quantize_coordinate(position[1], translate[1], scale[1])
      ]
      quantized.concat(position[2..-1]) if position.length > 2
      quantized
    end

    def quantize_coordinate(value, translate, scale)
      quantized = ((value - translate) / scale + 0.5).floor
      [[quantized, 0].max, @quantization - 1].min
    end

    def ring_invalid?(ring, ring_index, context)
      return false if ring.length >= 4 && !ring_area(ring).zero?
      diagnostic(:degenerate_ring, context.merge(ring: ring_index))
      true
    end

    def ring_area(ring)
      ring.each_cons(2).sum do |first, second|
        first[0] * second[1] - second[0] * first[1]
      end
    end

    def base_context(object, name: nil)
      context = {}
      context[:object] = name unless name.nil?
      context[:id] = object['id'] unless object['id'].nil?
      unless object['properties'].nil?
        context[:properties] = object['properties']
      end
      context
    end

    def child_context(child, parent_context, index)
      base_context(child).merge(parent_context).merge(
        geometry_path: (parent_context[:geometry_path] || []) + [index]
      )
    end

    def diagnostic(reason, context = {})
      @diagnostics << context.merge(reason: reason)
    end

    def build_topology(objects, transform)
      node_paths(objects)
      nodes, node_ids, edges, line_endpoints = build_graph(objects.values)
      @nodes = nodes
      @node_ids = node_ids
      replace_position_ids(objects, node_ids)
      arcs = build_arcs(nodes, node_ids, edges, line_endpoints)
      @arcs = arcs
      @arc_lengths = arcs.map(&:length)
      encoded = encode_arcs(arcs)
      edge_to_arc = register_arc_edges(arcs)
      exported = objects.transform_values do |object|
        clean_object(replace_paths(object, edge_to_arc))
      end
      {
        'type' => 'Topology',
        'objects' => exported,
        'arcs' => encoded,
        'transform' => transform
      }
    end

    def node_paths(objects)
      paths = []
      objects.each_value { |object| collect_paths(object, paths) }
      split = split_paths(paths)
      objects.each_value { |object| assign_paths(object, split) }
    end

    def collect_paths(object, paths)
      case object['type']
      when 'LineString' then paths << { positions: object['arcs'], open: true }
      when 'MultiLineString'
        object['arcs'].each { |line| paths << { positions: line, open: true } }
      when 'Polygon'
        object['arcs'].each { |ring| paths << { positions: ring, open: false } }
      when 'MultiPolygon'
        object['arcs'].each do |polygon|
          polygon.each { |ring| paths << { positions: ring, open: false } }
        end
      when 'GeometryCollection'
        object['geometries'].each { |child| collect_paths(child, paths) }
      end
    end

    def split_paths(paths)
      nodes = {}
      paths.each do |path|
        path[:positions].each { |position| nodes[position] = true }
      end
      nodes = nodes.keys
      paths.map do |path|
        split_path = []
        path[:positions].each_cons(2).with_index do |(start_node, end_node), edge_index|
          split_path << start_node if edge_index.zero?
          split_path.concat(nodes.select do |node|
            node != start_node && node != end_node &&
              point_on_segment?(node, start_node, end_node)
          end.sort_by { |node| segment_position(node, start_node, end_node) })
          split_path << end_node
        end
        split_path
      end
    end

    def point_on_segment?(node, start_node, end_node)
      cross = (node[1] - start_node[1]) * (end_node[0] - start_node[0]) -
              (node[0] - start_node[0]) * (end_node[1] - start_node[1])
      return false unless cross.zero?
      return false unless node.length == start_node.length &&
        node[2..-1] == start_node[2..-1] &&
        node[2..-1] == end_node[2..-1]
      between?(node[0], start_node[0], end_node[0]) &&
        between?(node[1], start_node[1], end_node[1])
    end

    def between?(value, start_value, end_value)
      [start_value, end_value].min <= value && value <= [start_value, end_value].max
    end

    def segment_position(node, start_node, end_node)
      delta_x = end_node[0] - start_node[0]
      delta_y = end_node[1] - start_node[1]
      (node[0] - start_node[0]) * delta_x +
        (node[1] - start_node[1]) * delta_y
    end

    def assign_paths(object, split)
      index = { value: 0 }
      replace_object_paths(object, index, split)
    end

    def replace_object_paths(object, index, split)
      case object['type']
      when 'LineString'
        object['arcs'] = split[index[:value]]
        index[:value] += 1
      when 'MultiLineString'
        object['arcs'] = object['arcs'].map do
          path = split[index[:value]]
          index[:value] += 1
          path
        end
      when 'Polygon'
        object['arcs'] = object['arcs'].map do
          path = split[index[:value]]
          index[:value] += 1
          path
        end
      when 'MultiPolygon'
        object['arcs'] = object['arcs'].map do |rings|
          rings.map do
            path = split[index[:value]]
            index[:value] += 1
            path
          end
        end
      when 'GeometryCollection'
        object['geometries'].each do |child|
          replace_object_paths(child, index, split)
        end
      end
    end

    def build_graph(objects)
      nodes = []
      node_ids = {}
      edges = {}
      line_endpoints = {}
      forced_nodes = {}
      paths = []
      objects.each { |object| collect_paths(object, paths) }
      paths.each do |path|
        path_nodes = path[:positions].map do |position|
          node_ids[position] ||= begin
            nodes << position
            nodes.length - 1
          end
        end
        line_endpoints[path_nodes.first] = true if path[:open]
        line_endpoints[path_nodes.last] = true if path[:open]
        path_nodes.each_cons(2) do |start_node, end_node|
          key = edge_key(start_node, end_node)
          edges[key] ||= [start_node, end_node]
        end
      end
      mark_noding_cuts(nodes, paths, node_ids, forced_nodes)
      [nodes, node_ids, edges, line_endpoints.merge(forced_nodes)]
    end

    def mark_noding_cuts(nodes, paths, node_ids, forced_nodes)
      vertices = nodes.each_with_object({}) { |node, set| set[node_ids[node]] = true }
      paths.each do |path|
        path[:positions].each_cons(2) do |start_node, end_node|
          start_id = node_ids[start_node]
          end_id = node_ids[end_node]
          nodes.each do |node|
            node_id = node_ids[node]
            next if node_id == start_id || node_id == end_id
            next unless vertices[node_id]
            forced_nodes[node_id] = true if point_on_segment?(
              node, start_node, end_node
            )
          end
        end
      end
    end

    def build_arcs(nodes, node_ids, edges, line_endpoints)
      neighbors = Hash.new { |hash, key| hash[key] = [] }
      incidence = Hash.new(0)
      edges.each_value do |start_node, end_node|
        next if start_node == end_node
        neighbors[start_node] << end_node
        neighbors[end_node] << start_node
        incidence[start_node] += 1
        incidence[end_node] += 1
      end
      cut = nodes.select do |node|
        index = node_ids[node]
        line_endpoints[index] || incidence[index] != 2
      end
      used = {}
      arcs = []
      cut.each do |node|
        node_id = node_ids[node]
        neighbors[node_id].sort.each do |neighbor|
          key = edge_key(node_id, neighbor)
          next if used[key]
          arcs << walk_arc(node_id, neighbor, used, neighbors, incidence,
                           cut_set(cut, node_ids))
        end
      end
      edges.each_value do |start_node, end_node|
        key = edge_key(start_node, end_node)
        next if used[key]
        arcs << walk_cycle(start_node, end_node, used, neighbors)
      end
      arcs
    end

    def cut_set(cut, node_ids)
      cut.each_with_object({}) { |node, set| set[node_ids[node]] = true }
    end

    def walk_arc(start_node, next_node, used, neighbors, incidence, cuts)
      positions = [start_node, next_node]
      used[edge_key(start_node, next_node)] = true
      current = next_node
      until cuts[current] || incidence[current] != 2
        following = neighbors[current].sort.find do |neighbor|
          !used[edge_key(current, neighbor)]
        end
        break if following.nil?
        break if cuts[following]
        used[edge_key(current, following)] = true
        positions << following
        current = following
      end
      positions
    end

    def walk_cycle(start_node, next_node, used, neighbors)
      positions = [start_node, next_node]
      used[edge_key(start_node, next_node)] = true
      previous = start_node
      current = next_node
      until current == start_node
        following = neighbors[current].sort.find do |node|
          node != previous && !used[edge_key(current, node)]
        end
        break if following.nil?
        used[edge_key(current, following)] = true
        positions << following
        previous = current
        current = following
      end
      positions
    end

    def edge_key(start_node, end_node)
      start_node < end_node ? [start_node, end_node] : [end_node, start_node]
    end

    def encode_arcs(arcs)
      arcs.map do |arc|
        positions = arc.map { |node_id| @nodes[node_id] }
        positions.map.with_index do |position, index|
          if index.zero?
            position
          else
            position.each_with_index.map do |coordinate, dimension|
              coordinate - positions[index - 1][dimension]
            end
          end
        end
      end
    end

    def register_arc_edges(arcs)
      mapping = {}
      arcs.each_with_index do |arc, arc_index|
        arc.each_cons(2).with_index do |(start_node, end_node), offset|
          mapping[[start_node, end_node]] = {
            arc: arc_index, offset: offset, reverse: false
          }
          mapping[[end_node, start_node]] = {
            arc: arc_index, offset: offset, reverse: true
          }
        end
      end
      mapping
    end

    def replace_position_ids(objects, node_ids)
      objects.each_value { |object| replace_object_position_ids(object, node_ids) }
    end

    def replace_object_position_ids(object, node_ids)
      type = object['type']
      if type == 'LineString'
        object['arcs'].map! { |position| node_ids[position] }
      elsif type == 'MultiLineString'
        object['arcs'].each do |line|
          line.map! { |position| node_ids[position] }
        end
      elsif type == 'Polygon'
        object['arcs'].each do |ring|
          ring.map! { |position| node_ids[position] }
        end
      elsif type == 'MultiPolygon'
        object['arcs'].each do |polygon|
          polygon.each do |ring|
            ring.map! { |position| node_ids[position] }
          end
        end
      elsif type == 'GeometryCollection'
        object['geometries'].each do |child|
          replace_object_position_ids(child, node_ids)
        end
      end
    end

    def clean_object(object)
      object.delete('_diagnostic')
      if object['type'] == 'GeometryCollection'
        object['geometries'].each { |child| clean_object(child) }
      end
      object
    end

    def replace_paths(object, mapping)
      case object['type']
      when 'Point', nil then object
      when 'LineString'
        object['arcs'] = line_arc_references(object['arcs'], mapping, true)
      when 'MultiLineString'
        object['arcs'] = object['arcs'].map do |line|
          line_arc_references(line, mapping, true)
        end
      when 'Polygon'
        object['arcs'] = object['arcs'].map do |ring|
          line_arc_references(ring, mapping, false)
        end
      when 'MultiPolygon'
        object['arcs'] = object['arcs'].map do |polygon|
          polygon.map { |ring| line_arc_references(ring, mapping, false) }
        end
      when 'GeometryCollection'
        object['geometries'] = object['geometries'].map do |child|
          replace_paths(child, mapping)
        end
      end
      object
    end

    def line_arc_references(positions, mapping, open_path)
      directed = positions.each_cons(2).map do |start_node, end_node|
        mapping.fetch([start_node, end_node])
      end
      if open_path
        return cyclic_references(directed) if positions.first == positions.last
        directed.map { |edge| arc_reference(edge) }
      else
        cyclic_references(directed)
      end
    end

    def cyclic_references(directed)
      return [] if directed.empty?
      boundary = directed.each_index.find do |index|
        next_index = (index + 1) % directed.length
        !consecutive_arc_edges?(directed[index], directed[next_index])
      end
      values = if boundary.nil?
                 directed
               elsif boundary == directed.length - 1
                 directed.rotate(boundary)
               else
                 directed.rotate(boundary + 1)
               end
      group_cyclic_edges(values)
    end

    def arc_reference(edge)
      edge[:reverse] ? ~edge[:arc] : edge[:arc]
    end

    def consecutive_arc_edges?(first, second)
      return false unless first[:arc] == second[:arc] &&
                            first[:reverse] == second[:reverse]
      edge_count = arc_edge_count(first[:arc])
      cyclic = cycle_arc?(first[:arc])
      if first[:reverse]
        expected = first[:offset] - 1
        expected = edge_count - 1 if cyclic && expected.negative?
      else
        expected = first[:offset] + 1
        expected = 0 if cyclic && expected == edge_count
      end
      expected == second[:offset]
    end

    def group_cyclic_edges(edges)
      groups = []
      edges.each do |edge|
        if groups.any? && append_to_cyclic_group?(groups.last, edge)
          groups.last << edge
        else
          groups << [edge]
        end
      end
      groups.map { |group| arc_reference(group.first) }
    end

    def append_to_cyclic_group?(group, edge)
      first = group.first
      return false unless consecutive_arc_edges?(group.last, edge)
      group.length < arc_edge_count(first[:arc])
    end

    def arc_edge_count(arc_index)
      @arc_lengths[arc_index] - 1
    end

    def cycle_arc?(arc_index)
      @arcs[arc_index].first == @arcs[arc_index].last
    end
  end
end
