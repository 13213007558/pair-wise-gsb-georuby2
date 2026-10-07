module GeoRuby
  module SimpleFeatures
    MeasurePoint = Struct.new(
      :point, :line_index, :segment_index, :parameter,
      :duplicate_vertex, :junction, keyword_init: true
    )

    MeasureSegment = Struct.new(
      :line, :line_index, :start_segment_index, :end_segment_index,
      :start_parameter, :end_parameter, :m_direction,
      :duplicate_vertex, :junction, keyword_init: true
    )

    MeasureSlice = Struct.new(
      :geometry, :kind, :line_index, :start_segment_index, :end_segment_index,
      :start_parameter, :end_parameter, :start_inclusive, :end_inclusive,
      :m_range, :m_direction, :duplicate_vertex, :junction, keyword_init: true
    )

    MissingMeasure = Struct.new(:line_index, :segment_index, :vertex_index, keyword_init: true)
    DuplicateVertex = Struct.new(:line_index, :segment_index, :vertex_index, :point, keyword_init: true)
    LineJunction = Struct.new(:line_indices, :segment_indices, :vertex_indices, :endpoints, :point, keyword_init: true)

    MeasureQueryResult = Struct.new(:geometry, :points, :segments, :issues, :junctions, keyword_init: true)
    MeasureRangeResult = Struct.new(:geometry, :slices, :issues, :junctions, keyword_init: true)

    class LinearReferencing
      Node = Struct.new(:point, :line_index, :segment_index, :parameter, :included_start, :included_end, :duplicate, keyword_init: true)
      Span = Struct.new(:line_index, :segment_index, :start_parameter, :end_parameter, keyword_init: true)

      def initialize(geometry)
        @geometry = geometry
        if geometry.is_a?(MultiLineString)
          @with_z = geometry.geometries.any?(&:with_z?)
          @with_m = geometry.geometries.any?(&:with_m?)
          geometry.with_z = @with_z unless geometry.with_z?
          geometry.with_m = @with_m unless geometry.with_m?
          @lines = geometry.geometries
        else
          @with_z = geometry.with_z?
          @with_m = geometry.with_m?
          @lines = [geometry]
        end
      end

      def locate_along(measure)
        prepare
        validate_measure(measure)

        points = []
        segments = []
        @lines.each_with_index do |line, line_index|
          next unless line.with_m? && line.length >= 2

          collect_constant_segments(line, line_index, measure, segments)
          collect_isolated_points(line, line_index, measure, points, segments)
        end

        annotate_junction(points, segments)
        MeasureQueryResult.new(
          geometry: query_geometry(points, segments),
          points: points,
          segments: segments,
          issues: @issues,
          junctions: @junctions
        )
      end

      def locate_between(start_measure, end_measure, include_start: true, include_end: true)
        prepare
        unless finite_measure?(start_measure) && finite_measure?(end_measure)
          raise ArgumentError, 'measure bounds must be finite numbers'
        end

        if start_measure == end_measure
          return zero_width_result(start_measure, include_start && include_end)
        end

        low, low_inclusive, high, high_inclusive =
          ordered_bounds(start_measure, end_measure, include_start, include_end)
        slices = []

        @lines.each_with_index do |line, line_index|
          next unless line.with_m? && line.length >= 2
          slices.concat(clip_line(line, line_index, low, high, low_inclusive, high_inclusive))
        end
        add_boundary_points(slices, low, high, low_inclusive, high_inclusive)

        annotate_slice_junctions(slices)
        MeasureRangeResult.new(
          geometry: range_geometry(slices),
          slices: slices,
          issues: @issues,
          junctions: @junctions
        )
      end

      def interpolate_point(start_point, end_point, parameter)
        point = Point.new(@geometry.srid, @with_z, @with_m)
        point.x = interpolate_value(start_point.x, end_point.x, parameter)
        point.y = interpolate_value(start_point.y, end_point.y, parameter)
        point.z = interpolate_value(start_point.z, end_point.z, parameter) if @with_z
        point.m = interpolate_value(start_point.m, end_point.m, parameter) if @with_m
        point
      end

      private

      def prepare
        return if @prepared

        @prepared = true
        @issues = []
        @junctions = []
        infer_line_flags
        endpoint_records = Hash.new { |hash, key| hash[key] = [] }

        @lines.each_with_index do |line, line_index|
          unless line && line.with_m?
            @issues << MissingMeasure.new(line_index: line_index)
            next
          end

          line.each_with_index do |point, vertex_index|
            unless point.with_m? && finite_measure?(point.m)
              segment_index = vertex_index.zero? ? nil : vertex_index - 1
              @issues << MissingMeasure.new(
                line_index: line_index,
                segment_index: segment_index,
                vertex_index: vertex_index
              )
            end
            record = endpoint_record(line, line_index, vertex_index)
            endpoint_records[position_key(point)] << record if record[:endpoint]
          end

          line.each_cons(2).with_index do |(start_point, end_point), segment_index|
            next unless same_position?(start_point, end_point)

            @issues << DuplicateVertex.new(
              line_index: line_index,
              segment_index: segment_index,
              vertex_index: segment_index,
              point: start_point
            )
          end
        end

        endpoint_records.each do |key, records|
          next if records.map { |record| record[:line_index] }.uniq.length < 2

          @junctions << LineJunction.new(
            line_indices: records.map { |record| record[:line_index] },
            segment_indices: records.map { |record| record[:segment_index] },
            vertex_indices: records.map { |record| record[:vertex_index] },
            endpoints: records.map { |record| record[:endpoint] },
            point: records.first[:point]
          )
          @junction_keys ||= {}
          @junction_keys[key] = true
        end
      end

      def endpoint_record(line, line_index, vertex_index)
        endpoint =
          if vertex_index.zero?
            :start
          elsif vertex_index == line.length - 1
            :end
          end
        segment_index =
          if endpoint == :end
            line.length - 2
          elsif endpoint == :start
            0
          end
        {
          line_index: line_index,
          vertex_index: vertex_index,
          segment_index: segment_index,
          endpoint: endpoint,
          point: line[vertex_index]
        }
      end

      def infer_line_flags
        @lines.each do |line|
          next unless line && line.length.positive? && line.first.is_a?(Point)

          line.with_z = true if !line.with_z? && line.points.all?(&:with_z?)
          line.with_m = true if !line.with_m? && line.points.all?(&:with_m?)
        end
        @with_z = @lines.any?(&:with_z?)
        @with_m = @lines.any?(&:with_m?)
        @geometry.with_z = @with_z unless @geometry.with_z?
        @geometry.with_m = @with_m unless @geometry.with_m?
      end

      def collect_constant_segments(line, line_index, measure, segments)
        run_nodes = nil
        line.each_cons(2).with_index do |(start_point, end_point), segment_index|
          constant = nonzero_edge?(start_point, end_point) &&
                     start_point.m == measure && end_point.m == measure
          if constant
            run_nodes ||= [make_node(start_point, line_index, segment_index, 0.0)]
            run_nodes << make_node(end_point, line_index, segment_index, 1.0)
          elsif run_nodes
            segments << measure_segment(run_nodes, line_index)
            run_nodes = nil
          end
        end
        segments << measure_segment(run_nodes, line_index) if run_nodes
      end

      def collect_isolated_points(line, line_index, measure, points, segments)
        platform_keys = segments.select { |segment| segment.line_index == line_index }
                             .flat_map { |segment| segment.line.points.map { |point| point_key(point) } }
        seen = {}

        line.each_with_index do |point, vertex_index|
          next unless valid_measure_point?(point) && point.m == measure
          next if platform_keys.include?(point_key(point))

          segment_index, parameter = vertex_location(line, vertex_index)
          add_measure_point(
            points, seen, point, line_index, segment_index, parameter,
            incident_duplicate?(line, vertex_index)
          )
        end

        line.each_cons(2).with_index do |(start_point, end_point), segment_index|
          next unless valid_measure_point?(start_point) && valid_measure_point?(end_point)
          next if same_position?(start_point, end_point)
          next if start_point.m == measure || end_point.m == measure
          next unless measure_between?(measure, start_point.m, end_point.m)

          parameter = (measure - start_point.m).to_f / (end_point.m - start_point.m)
          add_measure_point(
            points, seen,
            interpolate_point(start_point, end_point, parameter),
            line_index, segment_index, parameter, false
          )
        end
      end

      def vertex_location(line, vertex_index)
        if vertex_index == line.length - 1
          [line.length - 2, 1.0]
        else
          [vertex_index, 0.0]
        end
      end

      def incident_duplicate?(line, vertex_index)
        indices = []
        indices << vertex_index - 1 if vertex_index.positive?
        indices << vertex_index if vertex_index < line.length - 1
        indices.any? do |segment_index|
          same_position?(line[segment_index], line[segment_index + 1])
        end
      end

      def add_measure_point(points, seen, point, line_index, segment_index, parameter, duplicate)
        key = point_key(point)
        return if seen[key]

        seen[key] = true
        points << MeasurePoint.new(
          point: point,
          line_index: line_index,
          segment_index: segment_index,
          parameter: parameter,
          duplicate_vertex: duplicate,
          junction: false
        )
      end

      def measure_segment(nodes, line_index)
        MeasureSegment.new(
          line: LineString.from_points(nodes.map(&:point), @geometry.srid, @with_z, @with_m),
          line_index: line_index,
          start_segment_index: nodes.first.segment_index,
          end_segment_index: nodes.last.segment_index,
          start_parameter: nodes.first.parameter,
          end_parameter: nodes.last.parameter,
          m_direction: :constant,
          duplicate_vertex: nodes.any?(&:duplicate),
          junction: false
        )
      end

      def clip_line(line, line_index, low, high, low_inclusive, high_inclusive)
        slices = []
        current_nodes = []
        current_spans = []

        flush_line = lambda do
          return if current_nodes.empty?

          slices << build_slice(current_nodes, :line, line_index)
          current_nodes = []
          current_spans = []
        end

        flush_point = lambda do |point, segment_index, parameter, included|
          flush_line.call
          node = make_node(point, line_index, segment_index, parameter, included, included, true)
          slices << build_slice([node], :point, line_index)
        end

        line.each_cons(2).with_index do |(start_point, end_point), segment_index|
          unless valid_measure_point?(start_point) && valid_measure_point?(end_point)
            flush_line.call
            next
          end

          if same_position?(start_point, end_point)
            if included_measure?(start_point.m, low, high, low_inclusive, high_inclusive)
              flush_point.call(start_point, segment_index, 0.0,
                               included_bound?(start_point.m, low, high, low_inclusive, high_inclusive))
            else
              flush_line.call
            end
            next
          end

          spans = edge_spans(start_point, end_point, segment_index, line_index,
                            low, high, low_inclusive, high_inclusive)
          if spans.empty?
            flush_line.call
            next
          end

          spans.each_with_index do |span, span_index|
            start_node, end_node = *span
            previous_span = current_spans.last
            can_continue = previous_span &&
                           nodes_connect?(current_nodes.last, start_node) &&
                           current_nodes.last.included_end && start_node.included_start &&
                           same_edge_direction?(previous_span, start_node, end_node)
            if can_continue
              current_nodes << end_node
            else
              flush_line.call
              current_nodes = [start_node, end_node]
            end
            current_spans << [segment_index, span_index, start_node, end_node]
          end
        end

        flush_line.call
        slices
      end

      def edge_spans(start_point, end_point, segment_index, line_index,
                     low, high, low_inclusive, high_inclusive)
        if start_point.m == end_point.m
          return [] unless included_measure?(start_point.m, low, high, low_inclusive, high_inclusive)
        end

        parameters = { 0.0 => true, 1.0 => true }
        [low, high].each do |bound|
          if crossing_measure?(start_point.m, end_point.m, bound)
            parameters[solve_parameter(start_point.m, end_point.m, bound)] = true
          end
        end

        ordered_parameters = parameters.keys.sort
        ordered_parameters.each_cons(2).map do |first_parameter, second_parameter|
          mid_parameter = (first_parameter + second_parameter) / 2.0
          mid_measure = interpolate_value(start_point.m, end_point.m, mid_parameter)
          next unless mid_measure.between?(low, high)

          first_included = parameter_bound_included?(
            first_parameter, start_point.m, end_point.m, low, high, low_inclusive, high_inclusive
          )
          second_included = parameter_bound_included?(
            second_parameter, start_point.m, end_point.m, low, high, low_inclusive, high_inclusive
          )
          [
            make_node(
              interpolate_point(start_point, end_point, first_parameter),
              line_index, segment_index, first_parameter, first_included, first_included
            ),
            make_node(
              interpolate_point(start_point, end_point, second_parameter),
              line_index, segment_index, second_parameter, second_included, second_included
            )
          ]
        end.compact
      end

      def build_slice(nodes, kind, line_index)
        geometry =
          if kind == :point
            nodes.first.point
          else
            LineString.from_points(nodes.map(&:point), @geometry.srid, @with_z, @with_m)
          end
        first_node = nodes.first
        last_node = nodes.last
        MeasureSlice.new(
          geometry: geometry,
          kind: kind,
          line_index: line_index,
          start_segment_index: first_node.segment_index,
          end_segment_index: last_node.segment_index,
          start_parameter: first_node.parameter,
          end_parameter: last_node.parameter,
          start_inclusive: first_node.included_start,
          end_inclusive: last_node.included_end,
          m_range: kind == :point ? [geometry.m, geometry.m] : geometry.m_range,
          m_direction: slice_direction(geometry, kind),
          duplicate_vertex: nodes.any?(&:duplicate),
          junction: false
        )
      end

      def slice_direction(geometry, kind)
        return nil if kind == :point
        return :constant if geometry.length < 2

        edge_measures = geometry.points.each_cons(2).map do |start_point, end_point|
          [start_point.m, end_point.m]
        end
        if edge_measures.all? { |start_measure, end_measure| start_measure == end_measure }
          :constant
        elsif edge_measures.all? { |start_measure, end_measure| start_measure < end_measure }
          :ascending
        elsif edge_measures.all? { |start_measure, end_measure| start_measure > end_measure }
          :descending
        else
          :mixed
        end
      end

      def zero_width_result(measure, included)
        query_result = locate_along(measure)
        slices = []
        if included
          query_result.points.each do |point_hit|
            slices << MeasureSlice.new(
              geometry: point_hit.point,
              kind: :point,
              line_index: point_hit.line_index,
              start_segment_index: point_hit.segment_index,
              end_segment_index: point_hit.segment_index,
              start_parameter: point_hit.parameter,
              end_parameter: point_hit.parameter,
              start_inclusive: true,
              end_inclusive: true,
              m_range: [point_hit.point.m, point_hit.point.m],
              m_direction: nil,
              duplicate_vertex: point_hit.duplicate_vertex,
              junction: point_hit.junction
            )
          end
          query_result.segments.each do |segment|
            slices << MeasureSlice.new(
              geometry: segment.line,
              kind: :line,
              line_index: segment.line_index,
              start_segment_index: segment.start_segment_index,
              end_segment_index: segment.end_segment_index,
              start_parameter: segment.start_parameter,
              end_parameter: segment.end_parameter,
              start_inclusive: true,
              end_inclusive: true,
              m_range: segment.line.m_range,
              m_direction: segment.m_direction,
              duplicate_vertex: segment.duplicate_vertex,
              junction: segment.junction
            )
          end
        end
        MeasureRangeResult.new(
          geometry: range_geometry(slices),
          slices: slices,
          issues: query_result.issues,
          junctions: query_result.junctions
        )
      end

      def add_boundary_points(slices, low, high, low_inclusive, high_inclusive)
        bounds = []
        bounds << low if low_inclusive
        bounds << high if high_inclusive
        bounds.uniq.each do |bound|
          locate_along(bound).points.each do |point_hit|
            next if slices.any? do |slice|
              next false if slice.kind == :point

              slice.geometry.points.any? { |point| point_key(point) == point_key(point_hit.point) }
            end

            slices << MeasureSlice.new(
              geometry: point_hit.point,
              kind: :point,
              line_index: point_hit.line_index,
              start_segment_index: point_hit.segment_index,
              end_segment_index: point_hit.segment_index,
              start_parameter: point_hit.parameter,
              end_parameter: point_hit.parameter,
              start_inclusive: true,
              end_inclusive: true,
              m_range: [point_hit.point.m, point_hit.point.m],
              m_direction: nil,
              duplicate_vertex: point_hit.duplicate_vertex,
              junction: point_hit.junction
            )
          end
        end
      end

      def ordered_bounds(start_measure, end_measure, include_start, include_end)
        if start_measure < end_measure
          [start_measure, include_start, end_measure, include_end]
        else
          [end_measure, include_end, start_measure, include_start]
        end
      end

      def query_geometry(points, segments)
        point_geometry =
          if points.length == 1
            points.first.point
          elsif points.length > 1
            MultiPoint.from_points(points.map(&:point), @geometry.srid, @with_z, @with_m)
          end
        line_geometry = lines_geometry(segments.map(&:line))
        geometries = [point_geometry, line_geometry].compact
        return nil if geometries.empty?
        return geometries.first if geometries.length == 1

        GeometryCollection.from_geometries(geometries, @geometry.srid, @with_z, @with_m)
      end

      def range_geometry(slices)
        points = slices.select { |slice| slice.kind == :point }.map(&:geometry)
        lines = slices.select { |slice| slice.kind == :line }.map(&:geometry)
        geometries = []
        geometries << (points.length == 1 ? points.first : MultiPoint.from_points(points, @geometry.srid, @with_z, @with_m)) unless points.empty?
        geometries << lines_geometry(lines) unless lines.empty?
        return nil if geometries.empty?
        return geometries.first if geometries.length == 1

        GeometryCollection.from_geometries(geometries, @geometry.srid, @with_z, @with_m)
      end

      def lines_geometry(lines)
        return nil if lines.empty?
        return lines.first if @lines.length == 1 && lines.length == 1

        MultiLineString.from_line_strings(lines, @geometry.srid, @with_z, @with_m)
      end

      def annotate_junction(points, segments)
        return unless @junction_keys

        points.each do |point|
          point.junction = true if @junction_keys[position_key(point.point)]
        end
        segments.each do |segment|
          segment.junction = segment.line.points.any? do |point|
            @junction_keys[position_key(point)]
          end
        end
      end

      def annotate_slice_junctions(slices)
        return unless @junction_keys

        slices.each do |slice|
          points = slice.kind == :point ? [slice.geometry] : slice.geometry.points
          slice.junction = points.any? { |point| @junction_keys[position_key(point)] }
        end
      end

      def make_node(point, line_index, segment_index, parameter, included_start = true, included_end = true, duplicate = false)
        Node.new(
          point: point,
          line_index: line_index,
          segment_index: segment_index,
          parameter: parameter,
          included_start: included_start,
          included_end: included_end,
          duplicate: duplicate
        )
      end

      def nodes_connect?(first_node, second_node)
        first_node.line_index == second_node.line_index &&
          first_node.segment_index + 1 == second_node.segment_index &&
          first_node.parameter == 1.0 && second_node.parameter == 0.0 &&
          same_position?(first_node.point, second_node.point)
      end

      def same_edge_direction?(previous_span, start_node, end_node)
        previous_start_node = previous_span[2]
        previous_end_node = previous_span[3]
        previous_direction = measure_direction(previous_start_node.point.m, previous_end_node.point.m)
        current_direction = measure_direction(start_node.point.m, end_node.point.m)
        previous_direction == current_direction
      end

      def measure_direction(start_measure, end_measure)
        return :constant if start_measure == end_measure

        start_measure < end_measure ? :ascending : :descending
      end

      def valid_measures?(*measures)
        measures.all? { |measure| measure.is_a?(Numeric) && measure.finite? }
      end

      def valid_measure_point?(point)
        point.with_m? && finite_measure?(point.m)
      end

      def finite_measure?(measure)
        measure.is_a?(Numeric) && measure.finite?
      end

      def validate_measure(measure)
        raise ArgumentError, 'measure must be a finite number' unless finite_measure?(measure)
      end

      def nonzero_edge?(start_point, end_point)
        valid_measure_point?(start_point) && valid_measure_point?(end_point) && !same_position?(start_point, end_point)
      end

      def same_position?(first_point, second_point)
        if @with_z
          first_point.x == second_point.x &&
            first_point.y == second_point.y &&
            first_point.z == second_point.z
        else
          first_point.x == second_point.x && first_point.y == second_point.y
        end
      end

      def measure_between?(measure, first_measure, second_measure)
        if first_measure < second_measure
          measure > first_measure && measure < second_measure
        else
          measure > second_measure && measure < first_measure
        end
      end

      def crossing_measure?(first_measure, second_measure, measure)
        (first_measure < measure && second_measure > measure) ||
          (first_measure > measure && second_measure < measure)
      end

      def solve_parameter(start_measure, end_measure, measure)
        (measure - start_measure).to_f / (end_measure - start_measure)
      end

      def included_measure?(measure, low, high, low_inclusive, high_inclusive)
        (low_inclusive ? measure >= low : measure > low) &&
          (high_inclusive ? measure <= high : measure < high)
      end

      def included_bound?(measure, low, high, low_inclusive, high_inclusive)
        return low_inclusive if measure == low
        return high_inclusive if measure == high

        true
      end

      def parameter_bound_included?(parameter, start_measure, end_measure, low, high, low_inclusive, high_inclusive)
        measure =
          if parameter.zero?
            start_measure
          elsif parameter == 1.0
            end_measure
          else
            interpolate_value(start_measure, end_measure, parameter)
          end
        included_bound?(measure, low, high, low_inclusive, high_inclusive)
      end

      def interpolate_value(start_value, end_value, parameter)
        start_value + (end_value - start_value) * parameter
      end

      def position_key(point)
        if @with_z
          [point.x.to_f, point.y.to_f, point.z.to_f]
        else
          [point.x.to_f, point.y.to_f]
        end
      end

      def point_key(point)
        if @with_z
          [point.x.to_f, point.y.to_f, point.z.to_f, point.m.to_f]
        else
          [point.x.to_f, point.y.to_f, point.m.to_f]
        end
      end
    end
  end
end
