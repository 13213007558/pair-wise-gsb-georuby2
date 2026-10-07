module GeoRuby
  module SimpleFeatures
    # Continuous planar Fréchet distance for LineString geometries.
    class ContinuousFrechet
      NUMERIC_TOLERANCE = 1e-12
      PARAMETER_TOLERANCE = 1e-12

      Result = Struct.new(:distance, :error_bound, :lower_bound,
                          :upper_bound, :iterations)

      Decision = Struct.new(:within, :threshold, :numeric_tolerance,
                            :maximum_evidence_distance, :evidence)
      class Decision
        def within?
          within
        end
      end

      Node = Struct.new(:u, :v, :interval)

      Interval = Struct.new(:lo, :hi, :cell_i, :cell_j, :side, :source)
      class Interval
        def node_at(parameter)
          parameter = clamp_parameter(parameter)
          case side
          when :right
            Node.new(cell_i + 1, cell_j + parameter, self)
          when :top
            Node.new(cell_i + parameter, cell_j + 1, self)
          else
            Node.new(0, 0, self)
          end
        end

        def cover?(parameter)
          parameter >= lo - PARAMETER_TOLERANCE &&
            parameter <= hi + PARAMETER_TOLERANCE
        end

        private

        def clamp_parameter(parameter)
          return lo if parameter < lo && parameter >= lo - PARAMETER_TOLERANCE
          return hi if parameter > hi && parameter <= hi + PARAMETER_TOLERANCE
          parameter
        end
      end

      def initialize(first_line, second_line, tolerance = NUMERIC_TOLERANCE)
        @first_line = first_line
        @second_line = second_line
        @tolerance = tolerance
        validate!
        @first = coordinates(first_line)
        @second = coordinates(second_line)
      end

      def distance(iterations = 64)
        if point?(@first) || point?(@second)
          value = point_to_line_distance
          error = numeric_slack(value)
          return Result.new(value, error, value, value + error, 0)
        end

        lower = 0.0
        upper = max_vertex_distance
        upper *= 2.0 until decision(upper).within?

        iterations.times do
          candidate = (lower + upper) / 2.0
          if decision(candidate).within?
            upper = candidate
          else
            lower = candidate
          end
        end

        lower_error = numeric_slack(lower)
        upper_error = numeric_slack(upper)
        lower_bound = [0.0, lower - lower_error].max
        upper_bound = upper + upper_error
        value = (lower + upper) / 2.0
        error = [value - lower_bound, upper_bound - value].max
        Result.new(value, error, lower_bound, upper_bound, iterations)
      end

      def decision(threshold)
        raise ArgumentError, 'threshold must be non-negative' if threshold.negative?

        radius = threshold + numeric_slack(threshold)
        if point?(@first) || point?(@second)
          return point_decision(threshold, radius)
        end

        evidence = free_space_path(radius)
        return failed_decision(threshold) unless evidence

        maximum_distance = evidence.map { |step| step[:distance] }.max
        Decision.new(true, threshold, numeric_slack(threshold),
                      maximum_distance, evidence)
      end

      private

      def validate!
        unless @first_line.is_a?(LineString) && @second_line.is_a?(LineString)
          raise ArgumentError, 'both arguments must be LineString instances'
        end
        if @first_line.with_z? || @second_line.with_z?
          raise ArgumentError, 'continuous Fréchet distance is planar only'
        end
        if @first_line.empty? || @second_line.empty?
          raise ArgumentError, 'line strings must contain at least one point'
        end
        raise ArgumentError, 'tolerance must be non-negative' if @tolerance.negative?
      end

      def coordinates(line)
        line.points.map { |point| [point.x.to_f, point.y.to_f] }
      end

      def point?(line)
        line.length == 1
      end

      def point_to_line_distance
        if point?(@first)
          point = @first.first
          @second.map { |other| euclidian_distance(point, other) }.max
        else
          point = @second.first
          @first.map { |other| euclidian_distance(point, other) }.max
        end.to_f
      end

      def point_decision(threshold, radius)
        if point?(@first)
          point = @first.first
          valid = @second.all? { |other| euclidian_distance(point, other) <= radius }
          evidence = valid ? point_line_evidence(@first, @second, true) : nil
        else
          point = @second.first
          valid = @first.all? { |other| euclidian_distance(point, other) <= radius }
          evidence = valid ? point_line_evidence(@first, @second, false) : nil
        end

        maximum_distance = evidence ? evidence.map { |step| step[:distance] }.max : nil
        Decision.new(valid, threshold, numeric_slack(threshold),
                      maximum_distance, evidence)
      end

      def point_line_evidence(first, second, first_is_point)
        if first_is_point
          second.each_with_index.map do |_, index|
            parameter = normalize_parameter(index, second.length - 1)
            evidence_step(first.first, second[index], 0.0, parameter)
          end
        else
          first.each_with_index.map do |_, index|
            parameter = normalize_parameter(index, first.length - 1)
            evidence_step(first[index], second.first, parameter, 0.0)
          end
        end
      end

      def failed_decision(threshold)
        Decision.new(false, threshold, numeric_slack(threshold), nil, nil)
      end

      def max_vertex_distance
        @first.product(@second).map do |point, other|
          euclidian_distance(point, other)
        end.max
      end

      def numeric_slack(value)
        value.abs * @tolerance + @tolerance
      end

      def euclidian_distance(point, other)
        Math.sqrt((point[0] - other[0])**2 + (point[1] - other[1])**2)
      end

      def free_space_path(radius)
        # One cell per segment pair: a threshold decision is O(m*n) time and
        # space, and the distance routine performs a polynomial number of them.
        start_distance = euclidian_distance(@first.first, @second.first)
        return nil if start_distance > radius

        seed = Interval.new(0.0, 0.0, nil, nil, :seed, nil)
        segment_count_first = @first.length - 1
        segment_count_second = @second.length - 1
        right_intervals = Array.new(segment_count_first) do
          Array.new(segment_count_second)
        end
        top_intervals = Array.new(segment_count_first) do
          Array.new(segment_count_second)
        end
        target = nil

        segment_count_first.times do |i|
          segment_count_second.times do |j|
            left = if i.zero? && j.zero?
                     seed
                   elsif i.zero?
                     nil
                   else
                     right_intervals[i - 1][j]
                   end
            bottom = if i.zero? && j.zero?
                       seed
                     elsif j.zero?
                       nil
                     else
                       top_intervals[i][j - 1]
                     end

            a0 = @first[i]
            a1 = @first[i + 1]
            b0 = @second[j]
            b1 = @second[j + 1]

            left_interval = segment_interval(b0, b1, a0, radius)
            right_interval = segment_interval(b0, b1, a1, radius)
            bottom_interval = segment_interval(a0, a1, b0, radius)
            top_interval = segment_interval(a0, a1, b1, radius)

            left_source = reachable_source(left, left_interval)
            bottom_source = reachable_source(bottom, bottom_interval)

            right = reachable_right(right_interval, left_source, bottom_source, i, j)
            top = reachable_top(top_interval, left_source, bottom_source, i, j)

            right_intervals[i][j] = right
            top_intervals[i][j] = top
          end
        end

        last_right = right_intervals.last.last
        last_top = top_intervals.last.last
        if last_right&.cover?(1.0)
          target = last_right.node_at(1.0)
        elsif last_top&.cover?(1.0)
          target = last_top.node_at(1.0)
        end
        return nil unless target

        build_evidence(trace(target))
      end

      def reachable_source(interval, physical_interval)
        return nil unless interval && physical_interval

        parameter = [interval.lo, physical_interval[0]].max
        upper = [interval.hi, physical_interval[1]].min
        return nil if parameter > upper + PARAMETER_TOLERANCE

        parameter = upper if parameter > upper
        [parameter, interval.node_at(parameter)]
      end

      def reachable_right(physical_interval, left_source, bottom_source, i, j)
        return nil unless physical_interval

        if bottom_source
          _parameter, source = bottom_source
          Interval.new(physical_interval[0], physical_interval[1], i, j, :right, source)
        elsif left_source
          parameter, source = left_source
          lower = [physical_interval[0], parameter].max
          return nil if lower > physical_interval[1] + PARAMETER_TOLERANCE

          Interval.new(lower, physical_interval[1], i, j, :right, source)
        end
      end

      def reachable_top(physical_interval, left_source, bottom_source, i, j)
        return nil unless physical_interval

        if left_source
          _parameter, source = left_source
          Interval.new(physical_interval[0], physical_interval[1], i, j, :top, source)
        elsif bottom_source
          parameter, source = bottom_source
          lower = [physical_interval[0], parameter].max
          return nil if lower > physical_interval[1] + PARAMETER_TOLERANCE

          Interval.new(lower, physical_interval[1], i, j, :top, source)
        end
      end

      def segment_interval(start_point, end_point, fixed_point, radius)
        dx = end_point[0] - start_point[0]
        dy = end_point[1] - start_point[1]
        px = start_point[0] - fixed_point[0]
        py = start_point[1] - fixed_point[1]
        a = dx * dx + dy * dy
        b = 2.0 * (px * dx + py * dy)
        c = px * px + py * py - radius * radius

        if a.zero?
          return c <= 0 ? [0.0, 1.0] : nil
        end

        discriminant = b * b - 4.0 * a * c
        if discriminant < 0.0
          return c <= 0 ? [0.0, 1.0] : nil
        end

        root = Math.sqrt(discriminant)
        root_lower = (-b - root) / (2.0 * a)
        root_upper = (-b + root) / (2.0 * a)
        return nil if root_upper < 0.0 || root_lower > 1.0

        lower = [[root_lower, 0.0].max, 1.0].min
        upper = [[root_upper, 0.0].max, 1.0].min
        lower <= upper + PARAMETER_TOLERANCE ? [lower, upper] : nil
      end

      def trace(node, path = [])
        return path unless node

        trace(node.interval.source, path)
        path << node
        path
      end

      def build_evidence(nodes)
        evidence = nodes.map do |node|
          first_point = position(@first, node.u)
          second_point = position(@second, node.v)
          first_parameter = normalize_parameter(node.u, @first.length - 1)
          second_parameter = normalize_parameter(node.v, @second.length - 1)
          evidence_step(first_point, second_point,
                         first_parameter, second_parameter)
        end

        monotonic = evidence.each_cons(2).all? do |left_step, right_step|
          right_step[:u] >= left_step[:u] - PARAMETER_TOLERANCE &&
            right_step[:v] >= left_step[:v] - PARAMETER_TOLERANCE
        end
        monotonic ? evidence : nil
      end

      def normalize_parameter(parameter, segment_count)
        segment_count.zero? ? 0.0 : parameter / segment_count.to_f
      end

      def evidence_step(first_point, second_point, u, v)
        {
          u: u,
          v: v,
          first_point: first_point.dup,
          second_point: second_point.dup,
          distance: euclidian_distance(first_point, second_point)
        }
      end

      def position(line, parameter)
        return line.first if parameter <= 0
        return line.last if parameter >= line.length - 1

        segment = parameter.floor
        local = parameter - segment
        start_point = line[segment]
        end_point = line[segment + 1]
        [
          start_point[0] + (end_point[0] - start_point[0]) * local,
          start_point[1] + (end_point[1] - start_point[1]) * local
        ]
      end
    end
  end
end
