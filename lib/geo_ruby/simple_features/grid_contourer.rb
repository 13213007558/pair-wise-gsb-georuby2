module GeoRuby
  module SimpleFeatures
    class GridContour
      attr_reader :level, :line_strings, :linear_rings

      def initialize(level, line_strings, linear_rings)
        @level = level
        @line_strings = line_strings
        @linear_rings = linear_rings
      end

      def to_a
        line_strings + linear_rings
      end
    end

    class RegularGridContourer
      CellChain = Struct.new(:points, :start_angle, :end_angle)
      Occurrence = Struct.new(:chain_index, :end_index, :angle, :shared)

      attr_reader :values, :x_origin, :y_origin, :cell_width, :cell_height, :srid

      def initialize(values, x_origin = 0, y_origin = 0,
                     cell_width = 1, cell_height = 1,
                     srid = DEFAULT_SRID, missing_value: nil)
        @values = values
        @x_origin = x_origin
        @y_origin = y_origin
        @cell_width = cell_width
        @cell_height = cell_height
        @srid = srid
        @missing_value = missing_value

        raise ArgumentError, 'values must contain at least two rows' if row_count < 2
        raise ArgumentError, 'every row must contain at least two values' if col_count < 2
        raise ArgumentError, 'all rows must have the same length' unless values.all? { |row| row.length == col_count }
        raise ArgumentError, 'cell_width must be positive' unless cell_width > 0
        raise ArgumentError, 'cell_height must be positive' unless cell_height > 0
      end

      def contour(level)
        builder = LevelBuilder.new(self, level)
        builder.build
      end

      def contours(levels)
        levels.map { |level| contour(level) }
      end

      def row_count
        values.length
      end

      def col_count
        values.first.length
      end

      def missing?(value)
        value.nil? ||
          value.equal?(@missing_value) ||
          (value.is_a?(Float) && value.nan?)
      end

      class LevelBuilder
        attr_reader :point_cache

        CORNER_EDGES = {
          0 => [:bottom, :left],
          1 => [:bottom, :right],
          2 => [:right, :top],
          3 => [:left, :top]
        }.freeze

        def initialize(grid, level)
          @grid = grid
          @level = level
          @point_cache = {}
          @shareable_points = {}
          @chains = []
          @internal_pairs = {}
          @unshared_nodes = {}
          @level_edge_attachments = Hash.new { |hash, key| hash[key] = [] }
          @level_edge_keys = {}
        end

        def build
          @grid.row_count.times.to_a[0...-1].each do |row|
            @grid.col_count.times.to_a[0...-1].each do |col|
              add_cell(row, col)
            end
          end
          build_level_edges

          occurrences = collect_occurrences
          pairs = pair_occurrences(occurrences)
          pairs.merge!(@internal_pairs)
          line_points, ring_points = trace_geometries(occurrences, pairs)

          line_strings = line_points.map do |points|
            LineString.from_points(points, @grid.srid)
          end
          linear_rings = ring_points.map do |points|
            LinearRing.from_points(points, @grid.srid)
          end
          GridContour.new(@level, line_strings, linear_rings)
        end

        private

        def add_cell(row, col)
          positions = [
            [row, col],
            [row, col + 1],
            [row + 1, col + 1],
            [row + 1, col]
          ]
          cell_values = positions.map { |r, c| @grid.values[r][c] }
          return if cell_values.any? { |value| @grid.missing?(value) }

          @row = row
          @col = col
          @cell_values = cell_values
          @crossings = {}
          @crossing_parameters = {}

          zero_corners = cell_values.each_index.select { |i| cell_values[i] == @level }
          strict_crossings

          if zero_corners.empty?
            add_nonzero_cell
          elsif zero_corners.length == 1
            add_single_zero_corner_cell(zero_corners.first)
          elsif zero_corners.length == 2 && opposite?(zero_corners)
            add_opposite_zero_corners_cell(zero_corners)
          elsif zero_corners.length == 2
            add_adjacent_zero_corners_cell(zero_corners)
          end
        end

        def strict_crossings
          edges = {
            bottom: [0, 1],
            right: [1, 2],
            top: [3, 2],
            left: [0, 3]
          }
          edges.each do |name, (first, second)|
            first_value = @cell_values[first]
            second_value = @cell_values[second]
            next if first_value == @level || second_value == @level
            next if first_value > @level == (second_value > @level)

            parameter = (@level - first_value).to_f / (second_value - first_value)
            @crossings[name] = edge_point(name, parameter)
            @crossing_parameters[name] = parameter
          end
        end

        def add_nonzero_cell
          positives = @cell_values.each_index.select { |i| @cell_values[i] > @level }

          if positives.length == 1 || positives.length == 3
            isolated = positives.length == 1 ? positives : (0..3).to_a - positives
            add_chain([@crossings[CORNER_EDGES[isolated.first][0]],
                       @crossings[CORNER_EDGES[isolated.first][1]]])
          elsif positives.length == 2 && !opposite?(positives)
            add_chain(@crossings.values)
          elsif positives.length == 2
            add_ambiguous_cell(positives)
          end
        end

        def add_ambiguous_cell(positives)
          a, d, c, b = @cell_values
          determinant = a - b - d + c
          saddle_value = (a * c - b * d).to_f / determinant

          if saddle_value.zero?
            x_parameter = -(d - a).to_f / determinant
            y_parameter = -(b - a).to_f / determinant
            saddle = point([:saddle, @row, @col],
                           x_at(@col, x_parameter), y_at(@row, y_parameter))
            bottom_to_saddle = add_chain([edge_point(:bottom, x_parameter), saddle])
            saddle_to_top = add_chain([saddle, edge_point(:top, x_parameter)])
            left_to_saddle = add_chain([edge_point(:left, y_parameter), saddle])
            saddle_to_right = add_chain([saddle, edge_point(:right, y_parameter)])
            pair_internal_chains(bottom_to_saddle, saddle_to_top)
            pair_internal_chains(left_to_saddle, saddle_to_right)
            return
          end

          isolated = saddle_value > @level ? (0..3).to_a - positives : positives
          isolated.each do |corner|
            first_edge, second_edge = CORNER_EDGES[corner]
            add_chain([@crossings[first_edge], @crossings[second_edge]])
          end
        end

        def add_single_zero_corner_cell(corner)
          if @crossings.length == 1
            point = corner_point(corner)
            crossing = @crossings.values.first
            add_chain([point, crossing], corner_angle(corner), nil)
          elsif @crossings.length == 2
            add_chain(@crossings.values)
          end
        end

        def add_opposite_zero_corners_cell(corners)
          nonzero = (0..3).to_a - corners
          first_sign = @cell_values[nonzero.first] > @level
          second_sign = @cell_values[nonzero.last] > @level
          return if first_sign == second_sign

          add_chain([corner_point(corners.first), corner_point(corners.last)],
                    corner_angle(corners.first), corner_angle(corners.last))
        end

        def add_adjacent_zero_corners_cell(corners)
          opposite_corners = (0..3).to_a - corners
          return if opposite_corners.all? { |corner| @cell_values[corner] > @level } ||
                    opposite_corners.none? { |corner| @cell_values[corner] > @level }

          edge, opposite_edge = adjacent_zero_edge(corners.sort)
          parameter = attachment_parameter(edge, opposite_edge)
          return unless parameter

          attachment = edge_point(edge, parameter)
          @unshared_nodes[attachment] = true
          @level_edge_attachments[edge_key(edge)] << parameter
          add_chain([attachment, @crossings[opposite_edge]])
        end

        def adjacent_zero_edge(corners)
          case corners
          when [0, 1] then [:bottom, :top]
          when [1, 2] then [:right, :left]
          when [2, 3] then [:top, :bottom]
          when [0, 3] then [:left, :right]
          end
        end

        def attachment_parameter(edge, opposite_edge)
          parameter = @crossing_parameters[opposite_edge]
          return unless parameter
          edge == :top ? 1 - parameter : parameter
        end

        def add_chain(points, start_angle = nil, end_angle = nil)
          @chains << CellChain.new(points, start_angle, end_angle)
          @chains.length - 1
        end

        def pair_internal_chains(first_index, second_index)
          @internal_pairs[[first_index, 1]] = Occurrence.new(second_index, 0, nil, false)
          @internal_pairs[[second_index, 0]] = Occurrence.new(first_index, 1, nil, false)
        end

        def build_level_edges
          @grid.row_count.times do |row|
            (@grid.col_count - 1).times do |col|
              add_level_edge(:h, row, col) if active_level_edge?(:h, row, col)
            end
          end
          (@grid.row_count - 1).times do |row|
            @grid.col_count.times do |col|
              add_level_edge(:v, row, col) if active_level_edge?(:v, row, col)
            end
          end
        end

        def active_level_edge?(orientation, row, col)
          endpoint_values = if orientation == :h
                              [@grid.values[row][col], @grid.values[row][col + 1]]
                            else
                              [@grid.values[row][col], @grid.values[row + 1][col]]
                            end
          return false unless endpoint_values.all? { |value| value == @level }

          incident_cells(orientation, row, col).any? do |cell_row, cell_col|
            values = cell_values_at(cell_row, cell_col)
            values.none? { |value| @grid.missing?(value) } &&
              values.any? { |value| value != @level }
          end
        end

        def incident_cells(orientation, row, col)
          cells = if orientation == :h
                    [[row - 1, col], [row, col]]
                  else
                    [[row, col - 1], [row, col]]
                  end
          cells.select do |cell_row, cell_col|
            cell_row >= 0 && cell_row < @grid.row_count - 1 &&
              cell_col >= 0 && cell_col < @grid.col_count - 1
          end
        end

        def cell_values_at(row, col)
          [
            @grid.values[row][col],
            @grid.values[row][col + 1],
            @grid.values[row + 1][col + 1],
            @grid.values[row + 1][col]
          ]
        end

        def add_level_edge(orientation, row, col)
          key = [orientation, row, col]
          return if @level_edge_keys.key?(key)

          @level_edge_keys[key] = true
          points = [grid_edge_point(orientation, row, col, 0)]
          @level_edge_attachments[key].sort.each do |parameter|
            points << grid_edge_point(orientation, row, col, parameter)
          end
          points << grid_edge_point(orientation, row, col, 1)
          add_chain(points)
        end

        def opposite?(corners)
          (corners.first - corners.last).abs == 2
        end

        def edge_point(edge, parameter)
          grid_edge_point(*edge_key(edge), parameter)
        end

        def edge_key(edge)
          case edge
          when :bottom then [:h, @row, @col]
          when :top then [:h, @row + 1, @col]
          when :left then [:v, @row, @col]
          when :right then [:v, @row, @col + 1]
          end
        end

        def grid_edge_point(orientation, row, col, parameter)
          return grid_corner_point(row, col) if parameter.zero?
          if parameter == 1
            return grid_corner_point(row + (orientation == :v ? 1 : 0),
                                    col + (orientation == :h ? 1 : 0))
          end

          if orientation == :h
            point([orientation, row, col, parameter],
                  x_at(col, parameter), y_at(row))
          else
            point([orientation, row, col, parameter],
                  x_at(col), y_at(row, parameter))
          end
        end

        def corner_point(corner)
          row = @row + (corner < 2 ? 0 : 1)
          col = @col + (corner == 1 || corner == 2 ? 1 : 0)
          grid_corner_point(row, col)
        end

        def grid_corner_point(row, col)
          point([:corner, row, col], x_at(col), y_at(row))
        end

        def point(key, x, y)
          @point_cache[key] ||= begin
            point = Point.from_x_y(x, y, @grid.srid)
            @shareable_points[point] = true if %i[corner h v].include?(key.first)
            point
          end
        end

        def x_at(col, offset = 0)
          @grid.x_origin + (col + offset) * @grid.cell_width
        end

        def y_at(row, offset = 0)
          @grid.y_origin + (row + offset) * @grid.cell_height
        end

        def corner_angle(corner)
          x_neighbor = corner == 0 || corner == 3 ? 1 : 0
          y_neighbor = corner < 2 ? 3 : 2
          x_sign = corner == 0 || corner == 3 ? 1 : -1
          y_sign = corner < 2 ? 1 : -1

          x_derivative = x_sign * (@cell_values[x_neighbor] - @cell_values[corner])
          y_derivative = y_sign * (@cell_values[y_neighbor] - @cell_values[corner])
          tangent_x = y_derivative
          tangent_y = -x_derivative
          if tangent_x + tangent_y < 0
            tangent_x = -tangent_x
            tangent_y = -tangent_y
          end

          Math.atan2(y_sign * tangent_y, x_sign * tangent_x)
        end

        def collect_occurrences
          occurrences = Hash.new { |hash, key| hash[key] = [] }
          @chains.each_with_index do |chain, chain_index|
            start_angle = chain.start_angle || angle_between(chain.points[0], chain.points[1])
            end_angle = chain.end_angle || opposite_angle(chain.points[-1], chain.points[-2])
            occurrences[chain.points.first] << Occurrence.new(
              chain_index, 0, start_angle, shareable_node?(chain.points.first)
            )
            occurrences[chain.points.last] << Occurrence.new(
              chain_index, 1, end_angle, shareable_node?(chain.points.last)
            )
          end
          occurrences
        end

        def pair_occurrences(occurrences)
          pairs = {}
          occurrences.each do |point, node_occurrences|
            next unless @shareable_points.key?(point)
            remaining = node_occurrences.select(&:shared).dup
            while remaining.length > 1
              best = nil
              best_cost = nil
              remaining.combination(2).each do |first, second|
                next if first.chain_index == second.chain_index
                cost = opposite_cost(first.angle, second.angle)
                if best_cost.nil? || cost < best_cost
                  best = [first, second]
                  best_cost = cost
                end
              end
              first, second = best
              pairs[occurrence_key(first)] = second
              pairs[occurrence_key(second)] = first
              remaining.delete(first)
              remaining.delete(second)
            end
          end
          pairs
        end

        def opposite_cost(first_angle, second_angle)
          difference = ((first_angle - second_angle + Math::PI) % (2 * Math::PI)) - Math::PI
          (Math::PI - difference.abs).abs
        end

        def shareable_node?(point)
          @shareable_points.key?(point) && !@unshared_nodes.key?(point)
        end

        def trace_geometries(occurrences, pairs)
          used = {}
          open_lines = []
          closed_rings = []

          open_starts = occurrences.values.flatten.reject do |occurrence|
            pairs.key?(occurrence_key(occurrence))
          end

          open_starts.each do |start|
            next if used.key?(occurrence_key(start))
            open_lines << trace(start, pairs, used)
          end

          occurrences.values.flatten.each do |start|
            next if used.key?(occurrence_key(start))
            closed_rings << trace_closed(start, pairs, used)
          end

          [open_lines, closed_rings]
        end

        def trace(start, pairs, used)
          points = oriented_points(start)
          used[occurrence_key(start)] = true
          current = start

          while (next_occurrence = pairs[terminal_key(current)])
            used[terminal_key(current)] = true
            break if used.key?(occurrence_key(next_occurrence))
            append_oriented(points, next_occurrence)
            used[occurrence_key(next_occurrence)] = true
            current = next_occurrence
          end
          used[terminal_key(current)] = true
          points
        end

        def trace_closed(start, pairs, used)
          points = oriented_points(start)
          used[occurrence_key(start)] = true
          current = start

          while (next_occurrence = pairs[terminal_key(current)])
            used[terminal_key(current)] = true
            break if occurrence_key(next_occurrence) == occurrence_key(start)
            append_oriented(points, next_occurrence)
            used[occurrence_key(next_occurrence)] = true
            current = next_occurrence
          end
          points << points.first
        end

        def oriented_points(occurrence)
          points = @chains[occurrence.chain_index].points
          occurrence.end_index.zero? ? points.dup : points.reverse
        end

        def append_oriented(points, occurrence)
          oriented = oriented_points(occurrence)
          points.concat(oriented[1..-1])
        end

        def terminal_key(occurrence)
          [occurrence.chain_index, occurrence.end_index.zero? ? 1 : 0]
        end

        def occurrence_key(occurrence)
          [occurrence.chain_index, occurrence.end_index]
        end

        def angle_between(first, second)
          Math.atan2(second.y - first.y, second.x - first.x)
        end

        def opposite_angle(first, second)
          angle_between(first, second) + Math::PI
        end
      end
    end
  end
end
