require 'json'

module GeoRuby
  class TopoJSONDecoder
    def self.decode(topology, object_name: nil, quantized: false)
      new(topology).decode(object_name: object_name, quantized: quantized)
    end

    def self.coordinates(topology, object_name: nil, quantized: false)
      new(topology).coordinates(object_name: object_name, quantized: quantized)
    end

    def initialize(topology)
      @topology = topology
      @transform = topology['transform']
      @arcs = topology['arcs'].map { |arc| decode_arc_deltas(arc) }
    end

    def decode(object_name: nil, quantized: false)
      object = select_object(object_name)
      decode_object(object, quantized)
    end

    def coordinates(object_name: nil, quantized: false)
      object = select_object(object_name)
      object_coordinates(object, quantized)
    end

    private

    def select_object(name)
      objects = @topology['objects']
      name ||= objects.keys.first
      objects.fetch(name)
    end

    def decode_arc_deltas(arc)
      point = [0, 0]
      arc.map do |delta|
        point = [point.first + delta.first, point.last + delta.last]
        point
      end
    end

    def decode_object(object, quantized)
      return nil if object.nil?

      case object['type']
      when 'Point', 'MultiPoint'
        result = { 'type' => object['type'], 'coordinates' => object_coordinates(object, quantized) }
      when 'LineString', 'Polygon', 'MultiLineString', 'MultiPolygon'
        result = { 'type' => object['type'], 'coordinates' => object_coordinates(object, quantized) }
      when 'GeometryCollection'
        result = {
          'type' => 'GeometryCollection',
          'geometries' => object['geometries'].map { |member| decode_object(member, quantized) }
        }
      end
      copy_identity(object, result) if result
      result
    end

    def object_coordinates(object, quantized)
      case object['type']
      when 'Point'
        decode_point(object['coordinates'], quantized)
      when 'MultiPoint'
        object['coordinates'].map { |point| decode_point(point, quantized) }
      when 'LineString'
        decode_arc_chain(object['arcs'], quantized)
      when 'Polygon'
        object['arcs'].map { |ring| close_ring(decode_arc_chain(ring, quantized)) }
      when 'MultiLineString'
        object['arcs'].map { |line| decode_arc_chain(line, quantized) }
      when 'MultiPolygon'
        object['arcs'].map do |polygon|
          polygon.map { |ring| close_ring(decode_arc_chain(ring, quantized)) }
        end
      when 'GeometryCollection'
        object['geometries'].map { |member| object_coordinates(member, quantized) }
      end
    end

    def decode_point(point, quantized)
      return point if quantized || @transform.nil?

      [
        point.first * @transform['scale'].first + @transform['translate'].first,
        point.last * @transform['scale'].last + @transform['translate'].last
      ]
    end

    def decode_arc_chain(references, quantized)
      points = []
      Array(references).each do |reference|
        arc = referenced_arc(reference)
        arc.each_with_index do |point, index|
          next if !points.empty? && index.zero?
          points << decode_point(point, quantized)
        end
      end
      points
    end

    def referenced_arc(reference)
      if reference >= 0
        @arcs.fetch(reference)
      else
        @arcs.fetch(~reference).reverse
      end
    end

    def close_ring(points)
      points.empty? || points.first == points.last ? points : points + [points.first]
    end

    def copy_identity(source, result)
      result['id'] = source['id'] if source.key?('id')
      result['properties'] = source['properties'] if source.key?('properties')
    end
  end
end
