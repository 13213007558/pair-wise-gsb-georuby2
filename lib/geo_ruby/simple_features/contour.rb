require 'geo_ruby/simple_features/line_string'
require 'geo_ruby/simple_features/linear_ring'
require 'geo_ruby/simple_features/polygon'

module GeoRuby
  module SimpleFeatures
    # Builds iso-lines from a regular scalar grid. Each cell is interpreted
    # with bilinear interpolation. Saddle cells use an asymptotic decider
    # rather than always choosing the same diagonal.
    class Contour
      # Open contours and closed contours. Closed contours are returned as
      # LinearRing instances; #polygons wraps each ring in a Polygon.
      Result = Struct.new(:open_lines, :closed_rings) do
        alias_method :lines, :open_lines
        alias_method :rings, :closed_rings

        def polygons
          closed_rings.map do |ring|
            Polygon.from_linear_rings([ring])
          end
        end

        def to_a
          [open_lines, closed_rings]
        end
      end

      Arc = Struct.new(:start_id, :end_id, :points)

      BOTTOM = 0
      RIGHT = 1
      TOP = 2
      LEFT = 3

      CORNERS = [
        [0, 0],
        [1, 0],
        [1, 1],
        [0, 1]
      ].freeze

      def self.from_grid(grid, level:, origin_x: 0.0, origin_y: 0.0,
                         cell_width: 1.0, cell_height: 1.0,
                         srid: DEFAULT_SRID, nodata: nil,
                         curve_segments: 16)
        new(grid, level, origin_x, origin_y, cell_width, cell_height,
            srid, nodata, curve_segments).build
      end

      def initialize(grid, level, origin_x, origin_y, cell_width, cell_height,
                     srid, nodata, curve_segments)
        @grid = grid
        @level = level
        @origin_x = origin_x.to_f
        @origin_y = origin_y.to_f
        @cell_width = cell_width.to_f
        @cell_height = cell_height.to_f
        @srid = srid
        @nodata = nodata
        @curve_segments = curve_segments

        @points = {}
        @arcs = []
        @ports = Hash.new { |hash, key| hash[key] = [] }
        @connections = []
      end

      def build
        build_arcs
        pair_ports
        trace_geometries
        Result.new(@open_lines, @closed_rings)
      end

      private

      def build_arcs
        @grid.each_with_index do |row, cell_row|
          next if row.nil?

          row.each_index do |cell_column|
            build_cell_arcs(cell_row, cell_column)
          end
        end

        deduplicate_coincident_arcs
        register_ports
      end

      def build_cell_arcs(cell_row, cell_column)
        values = CORNERS.map do |corner_column, corner_row|
          value_at(cell_row + corner_row, cell_column + corner_column)
        end
        return if values.any? { |value| missing?(value) }

        bits = values.each_with_index.sum do |value, index|
          value >= @level ? 1 << index : 0
        end
        return if bits.zero? || bits == 15

        edges = edge_intersections(cell_row, cell_column, values)
        exact_saddle = [5, 10].include?(bits) && saddle_center(values) == @level
        saddle_id = cell_center_id(cell_row, cell_column) if exact_saddle
        pairings_for(bits, values).each do |pair|
          first_edge, first_index, second_edge, second_index = pair
          interior_points = bilinear_path_points(
            cell_row, cell_column, values, pair, exact_saddle
          )
          add_edge_arc(edges, first_edge, first_index,
                       second_edge, second_index, saddle_id, interior_points)
        end
      end

      def pairings_for(bits, values)
        case bits
        when 1 then [[BOTTOM, 0, LEFT, 0]]
        when 2 then [[BOTTOM, 1, RIGHT, 0]]
        when 3 then [[LEFT, 0, RIGHT, 0]]
        when 4 then [[RIGHT, 1, TOP, 1]]
        when 6 then [[BOTTOM, 1, TOP, 1]]
        when 7 then [[LEFT, 0, TOP, 1]]
        when 8 then [[LEFT, 1, TOP, 0]]
        when 9 then [[BOTTOM, 0, TOP, 0]]
        when 11 then [[RIGHT, 1, TOP, 0]]
        when 12 then [[LEFT, 0, RIGHT, 1]]
        when 13 then [[BOTTOM, 0, RIGHT, 0]]
        when 14 then [[BOTTOM, 1, LEFT, 1]]
        when 5
          center = saddle_center(values)
          if center > @level
            [[BOTTOM, 1, RIGHT, 0], [LEFT, 1, TOP, 0]]
          elsif center < @level
            [[BOTTOM, 0, LEFT, 0], [RIGHT, 1, TOP, 1]]
          else
            [[BOTTOM, 0, TOP, 0], [LEFT, 0, RIGHT, 0]]
          end
        when 10
          center = saddle_center(values)
          if center > @level
            [[BOTTOM, 0, LEFT, 0], [RIGHT, 1, TOP, 1]]
          elsif center < @level
            [[BOTTOM, 1, RIGHT, 0], [LEFT, 1, TOP, 0]]
          else
            [[BOTTOM, 0, TOP, 0], [LEFT, 0, RIGHT, 0]]
          end
        end
      end

      def saddle_center(values)
        (values[0] + values[1] + values[2] + values[3]) / 4.0
      end

      def edge_intersections(cell_row, cell_column, values)
        [
          intersection([:horizontal, cell_row, cell_column],
                       values[0], values[1],
                       [:node, cell_row, cell_column],
                       [:node, cell_row, cell_column + 1]),
          intersection([:vertical, cell_row, cell_column + 1],
                       values[1], values[2],
                       [:node, cell_row, cell_column + 1],
                       [:node, cell_row + 1, cell_column + 1]),
          intersection([:horizontal, cell_row + 1, cell_column],
                       values[3], values[2],
                       [:node, cell_row + 1, cell_column],
                       [:node, cell_row + 1, cell_column + 1]),
          intersection([:vertical, cell_row, cell_column],
                       values[0], values[3],
                       [:node, cell_row, cell_column],
                       [:node, cell_row + 1, cell_column])
        ]
      end

      def intersection(key, first_value, second_value, first_node, second_node)
        difference = second_value.to_f - first_value.to_f
        if difference.zero?
          return nil unless first_value == @level

          return [:flat_edge, first_node, second_node]
        end

        parameter = (@level - first_value) / difference
        return nil unless parameter >= 0.0 && parameter <= 1.0

        if parameter.zero? || parameter == 1.0
          parameter.zero? ? first_node : second_node
        elsif parameter > 0.0 && parameter < 1.0
          if key[0] == :horizontal
            x = @origin_x + (grid_column(key) + parameter) * @cell_width
            y = @origin_y + grid_row(key) * @cell_height
          else
            x = @origin_x + grid_column(key) * @cell_width
            y = @origin_y + (grid_row(key) + parameter) * @cell_height
          end
          point_id(key, x, y)
        else
          nil
        end
      end

      def grid_row(key)
        key[1]
      end

      def grid_column(key)
        key[2]
      end

      def point_id(key, x, y)
        id = [:point, key, x, y]
        @points[id] ||= Point.from_x_y(x, y, @srid)
        id
      end

      def cell_center_id(cell_row, cell_column)
        [:cell_center, cell_row, cell_column]
      end

      def cell_center_point(id)
        @points[id] ||= Point.from_x_y(
          @origin_x + (id[2] + 0.5) * @cell_width,
          @origin_y + (id[1] + 0.5) * @cell_height,
          @srid
        )
      end

      def bilinear_path_points(cell_row, cell_column, values, pair, exact_saddle)
        return [] if exact_saddle || @curve_segments <= 1

        first_edge, _first_index, second_edge, _second_index = pair
        return [] if first_edge == second_edge

        horizontal_edge = [first_edge, second_edge].find do |edge|
          [BOTTOM, TOP].include?(edge)
        end
        vertical_edge = [first_edge, second_edge].find do |edge|
          [LEFT, RIGHT].include?(edge)
        end

        if horizontal_edge && vertical_edge
          points = mixed_edge_points(values, horizontal_edge, vertical_edge)
        elsif vertical_edge
          points = vertical_edge_points(values, pair)
        else
          points = horizontal_edge_points(values, pair)
        end

        points.reverse! if first_edge == vertical_edge
        points.map do |u, v|
          x = @origin_x + (cell_column + u) * @cell_width
          y = @origin_y + (cell_row + v) * @cell_height
          Point.from_x_y(x, y, @srid)
        end
      end

      def mixed_edge_points(values, horizontal_edge, vertical_edge)
        start_u = edge_parameter(horizontal_edge, values)
        start_v = horizontal_edge == TOP ? 1.0 : 0.0
        end_u = vertical_edge == RIGHT ? 1.0 : 0.0
        end_v = edge_parameter(vertical_edge, values)

        (1...@curve_segments).map do |index|
          fraction = index.to_f / @curve_segments
          u = start_u + (end_u - start_u) * fraction
          v = bilinear_v(values, u, start_v, end_v)
          [u, v]
        end
      end

      def vertical_edge_points(values, pair)
        first_edge, _first_index, second_edge, _second_index = pair
        first_parameter = edge_parameter(first_edge, values)
        second_parameter = edge_parameter(second_edge, values)

        (1...@curve_segments).map do |index|
          fraction = index.to_f / @curve_segments
          v = first_parameter + (second_parameter - first_parameter) * fraction
          [bilinear_u(values, v, first_parameter, second_parameter), v]
        end
      end

      def horizontal_edge_points(values, pair)
        first_edge, _first_index, second_edge, _second_index = pair
        first_parameter = edge_parameter(first_edge, values)
        second_parameter = edge_parameter(second_edge, values)

        (1...@curve_segments).map do |index|
          fraction = index.to_f / @curve_segments
          u = first_parameter + (second_parameter - first_parameter) * fraction
          v = bilinear_horizontal_v(values, u) do
            first_parameter + (second_parameter - first_parameter) * fraction
          end
          [u, v]
        end
      end

      def edge_parameter(edge, values)
        case edge
        when BOTTOM then parameter_between(values[0], values[1])
        when RIGHT then parameter_between(values[1], values[2])
        when TOP
          parameter_between(values[3], values[2])
        when LEFT then parameter_between(values[0], values[3])
        end
      end

      def parameter_between(first_value, second_value)
        difference = second_value.to_f - first_value.to_f
        return 0.5 if difference.zero?

        parameter = (@level - first_value) / difference
        [[parameter, 0.0].max, 1.0].min
      end

      def bilinear_v(values, u, start_v, end_v)
        a = values[0]
        b = values[1]
        c = values[2]
        d = values[3]
        denominator = d - a + (a + c - b - d) * u

        if denominator.zero?
          return start_v + (end_v - start_v) * u
        end

        numerator = @level - a - (b - a) * u
        linear_v = numerator / denominator
        [[linear_v, 0.0].max, 1.0].min
      end

      def bilinear_u(values, v, start_u, end_u)
        a = values[0]
        b = values[1]
        c = values[2]
        d = values[3]
        denominator = b - a + (a + d - b - c) * v
        return start_u + (end_u - start_u) * v if denominator.zero?

        numerator = @level - a - (d - a) * v
        linear_u = numerator / denominator
        [[linear_u, 0.0].max, 1.0].min
      end

      def bilinear_horizontal_v(values, u)
        a = values[0]
        b = values[1]
        c = values[2]
        d = values[3]
        denominator = d - a + (a + c - b - d) * u
        return yield if denominator.zero?

        numerator = @level - a - (b - a) * u
        linear_v = numerator / denominator
        [[linear_v, 0.0].max, 1.0].min
      end

      def add_edge_arc(edges, first_edge, first_index,
                       second_edge, second_index,
                       interior_id = nil, interior_points = [])
        first_edge_points = edges[first_edge]
        second_edge_points = edges[second_edge]
        first_id = edge_point(first_edge_points, first_index)
        second_id = edge_point(second_edge_points, second_index)
        return if first_id.is_a?(Array) && first_id.first.nil?
        return if second_id.is_a?(Array) && second_id.first.nil?
        return if first_id.nil? || second_id.nil?

        add_arc(first_id, second_id, interior_id, interior_points)
      end

      def edge_point(edge, corner_index)
        return nil if edge.nil?
        return edge[corner_index + 1] if edge.first == :flat_edge

        edge
      end

      def add_arc(first_id, second_id, interior_id = nil, interior_points = [])
        return if first_id.nil? || second_id.nil?
        return if first_id == second_id

        points = [point_for(first_id)]
        points.concat(interior_points)
        points << point_for(interior_id) if interior_id
        points << point_for(second_id)
        points = compact_points(points)
        return if points.length < 2
        if interior_id
          first_points = [points.first, point_for(interior_id)]
          second_points = [point_for(interior_id), points.last]
          first_arc = Arc.new(first_id, interior_id, first_points)
          second_arc = Arc.new(interior_id, second_id, second_points)
          @arcs << first_arc << second_arc
        else
          @arcs << Arc.new(first_id, second_id, points)
        end
      end

      def compact_points(points)
        points.chunk do |point|
          [point.x.round(12), point.y.round(12)]
        end.map(&:last).map(&:first)
      end

      def point_for(id)
        return @points[id] if id[0] == :point
        return cell_center_point(id) if id[0] == :cell_center

        @points[id] ||= Point.from_x_y(
          @origin_x + id[2] * @cell_width,
          @origin_y + id[1] * @cell_height,
          @srid
        )
      end

      def deduplicate_coincident_arcs
        seen = {}
        @arcs.select! do |arc|
          key = [arc.start_id, arc.end_id].minmax
          !seen.key?(key) && seen[key] = true
        end
      end

      def register_ports
        @arcs.each_with_index do |arc, arc_index|
          register_port(arc, arc_index, 0, angle_from(arc.points[0], arc.points[1]))
          register_port(arc, arc_index, 1, angle_from(arc.points[1], arc.points[0]))
        end
      end

      def register_port(arc, arc_index, endpoint_index, angle)
        node_id = endpoint_index.zero? ? arc.start_id : arc.end_id
        @ports[node_id] << [arc_index, endpoint_index, angle]
      end

      def angle_from(first_point, second_point)
        Math.atan2(second_point.y - first_point.y,
                   second_point.x - first_point.x)
      end

      def pair_ports
        @connections = Array.new(@arcs.size) { [nil, nil] }

        @ports.each_value do |ports|
          next unless ports.size == 2 || ports.size == 4

          ordered = ports.sort_by { |port| [port[2], port[0], port[1]] }
          pairs = ports.size == 2 ? [[ordered[0], ordered[1]]] :
            [[ordered[0], ordered[2]], [ordered[1], ordered[3]]]
          pairs.each { |first, second| connect(first, second) }
        end
      end

      def connect(first, second)
        @connections[first[0]][first[1]] = [second[0], second[1]]
        @connections[second[0]][second[1]] = [first[0], first[1]]
      end

      def trace_geometries
        @open_lines = []
        @closed_rings = []
        consumed = Array.new(@arcs.size) { [false, false] }

        @arcs.each_with_index do |_arc, arc_index|
          [0, 1].each do |endpoint_index|
            next if consumed[arc_index][endpoint_index]
            next if @connections[arc_index][endpoint_index]

            @open_lines << LineString.from_points(
              trace_path(arc_index, endpoint_index, consumed), @srid
            )
          end
        end

        @arcs.each_with_index do |_arc, arc_index|
          next if consumed[arc_index][0] && consumed[arc_index][1]

          first_unconsumed = consumed[arc_index][0] ? 1 : 0
          points = trace_cycle(arc_index, first_unconsumed, consumed)
          next if points.length < 3

          @closed_rings << LinearRing.from_points(points, @srid)
        end
      end

      def trace_path(first_arc, first_endpoint, consumed)
        points = [@arcs[first_arc].points[first_endpoint]]
        arc_index = first_arc
        endpoint_index = first_endpoint

        loop do
          other_endpoint = 1 - endpoint_index
          points.concat(segment_points(arc_index, endpoint_index, points.last))
          consumed[arc_index][endpoint_index] = true
          consumed[arc_index][other_endpoint] = true

          next_port = @connections[arc_index][other_endpoint]
          break unless next_port

          arc_index, endpoint_index = next_port
        end

        points
      end

      def trace_cycle(first_arc, first_endpoint, consumed)
        points = [@arcs[first_arc].points[first_endpoint]]
        arc_index = first_arc
        endpoint_index = first_endpoint
        start_port = [first_arc, first_endpoint]

        loop do
          other_endpoint = 1 - endpoint_index
          consumed[arc_index][endpoint_index] = true
          consumed[arc_index][other_endpoint] = true

          next_port = @connections[arc_index][other_endpoint]
          points.concat(segment_points(arc_index, endpoint_index, points.last))
          break unless next_port
          break if next_port == start_port

          arc_index, endpoint_index = next_port
        end

        points
      end

      def segment_points(arc_index, endpoint_index, previous_point)
        arc = @arcs[arc_index]
        segment = if endpoint_index.zero?
                    arc.points[1..-1]
                  else
                    arc.points[1..-2].reverse + [arc.points[0]]
                  end
        segment.reject do |point|
          point.x == previous_point.x && point.y == previous_point.y
        end
      end

      def value_at(row, column)
        row_values = @grid[row]
        return nil if row_values.nil? || column < 0 || column >= row_values.length

        row_values[column]
      end

      def missing?(value)
        value.nil? || value.equal?(@nodata) || (value.respond_to?(:nan?) && value.nan?)
      end
    end
  end
end
