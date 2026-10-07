require 'geo_ruby/simple_features/geometry'

module GeoRuby
  module SimpleFeatures
    # Linear referencing on measured (M-aware) lines.
    #
    # The M ordinate is treated as a milepost: it is not necessarily
    # related to the coordinate length of the line and is not required
    # to be monotonic. All interpolation is done segment by segment,
    # against the M values of the two endpoints of each segment - never
    # against the whole line length. When the line walks back (the same
    # measure occurs several times) every occurrence is reported, and
    # when a whole segment keeps a constant M (a plateau) the segment
    # itself belongs to that measure.
    #
    # Included in LineString (and therefore LinearRing) and
    # MultiLineString.
    module LinearReferencing
      # Returns every location of the line where the interpolated M
      # value equals +m+.
      #
      # Isolated crossings come back as points, segments whose M is
      # constantly +m+ come back as lines. The result is a Point,
      # MultiPoint, LineString, MultiLineString or GeometryCollection
      # (when points and lines are mixed), or nil when the measure is
      # not found on the line.
      def locate_along(m)
        locate_between(m, m)
      end

      # Clips the line to the measures between +m_start+ and +m_end+
      # (in any order) and returns every sub-line inside the interval.
      #
      # +include_start+ and +include_end+ control whether the interval
      # endpoints are closed (true, the default) or open (false). The
      # distinction matters for degenerate results: an isolated point
      # sitting exactly on an open endpoint is left out, as is a whole
      # constant-M segment whose measure equals an open endpoint.
      #
      # Consecutive pieces sharing a junction point are merged, so the
      # result is a LineString when the clipped parts are connected, a
      # MultiLineString when they are not, a Point/MultiPoint when only
      # isolated locations match, a GeometryCollection for a mix, or
      # nil when nothing matches.
      #
      # Z (when present) is interpolated with the same parameter as the
      # X/Y coordinates. SRID, Z/M flags and the original direction of
      # travel are preserved in the output.
      def locate_between(m_start, m_end, include_start: true, include_end: true)
        unless with_m
          raise ArgumentError, 'geometry has no M ordinate to locate on'
        end
        lo, hi = m_start <= m_end ? [m_start, m_end] : [m_end, m_start]
        builder = MeasureResultBuilder.new(srid, with_z, with_m)
        each_measured_segment do |p1, p2|
          clip_segment_to_measures(p1, p2, lo, hi, include_start, include_end, builder)
        end
        builder.result
      end

      private

      # Yields each pair of consecutive vertices of the line (or of each
      # part of a multi line) in the original direction of travel.
      def each_measured_segment
        parts = is_a?(MultiLineString) ? geometries : [self]
        parts.each do |part|
          part.points.each_cons(2) { |p1, p2| yield p1, p2 }
        end
      end

      def clip_segment_to_measures(p1, p2, lo, hi, inc_lo, inc_hi, builder)
        m1 = p1.m
        m2 = p2.m
        if m1.nil? || m2.nil?
          # missing measure: the line cannot be interpolated here
          builder.break!
          return
        end
        m1 = m1.to_f
        m2 = m2.to_f
        zero_length = p1.x == p2.x && p1.y == p2.y
        if m1 == m2
          clip_plateau(p1, p2, m1, lo, hi, inc_lo, inc_hi, zero_length, builder)
        elsif zero_length
          # duplicate vertices spanning a measure jump: a single location
          clip_zero_length_spike(p1, m1, m2, lo, hi, inc_lo, inc_hi, builder)
        else
          clip_sloped_segment(p1, p2, m1, m2, lo, hi, inc_lo, inc_hi, builder)
        end
      end

      # A segment whose M is constant belongs entirely to one measure.
      def clip_plateau(p1, p2, m, lo, hi, inc_lo, inc_hi, zero_length, builder)
        unless measure_covered?(m, lo, hi, inc_lo, inc_hi)
          builder.break!
          return
        end
        if zero_length
          builder.add_point(p1)
        else
          builder.add_piece(p1, p2)
        end
      end

      def clip_zero_length_spike(p1, m1, m2, lo, hi, inc_lo, inc_hi, builder)
        m_lo, m_hi = m1 < m2 ? [m1, m2] : [m2, m1]
        o_lo = lo > m_lo ? lo : m_lo
        o_hi = hi < m_hi ? hi : m_hi
        if o_lo < o_hi || (o_lo == o_hi && boundary_closed?(o_lo, lo, hi, inc_lo, inc_hi))
          point = Point.new(srid, with_z, with_m)
          point.set_x_y(p1.x, p1.y)
          point.z = p1.z if with_z
          point.m = o_lo
          builder.add_point(point)
        else
          builder.break!
        end
      end

      # A segment whose M varies: interpolate the overlap of its measure
      # range with the interval, using the segment M values only.
      def clip_sloped_segment(p1, p2, m1, m2, lo, hi, inc_lo, inc_hi, builder)
        t_lo = (lo - m1) / (m2 - m1)
        t_hi = (hi - m1) / (m2 - m1)
        t_a, t_b = t_lo < t_hi ? [t_lo, t_hi] : [t_hi, t_lo]
        t_a = 0.0 if t_a < 0.0
        t_b = 1.0 if t_b > 1.0
        if t_a < t_b
          builder.add_piece(point_at_t(p1, p2, t_a, m1, m2, lo, hi),
                            point_at_t(p1, p2, t_b, m1, m2, lo, hi))
        else
          m_val = m1 + t_a * (m2 - m1)
          if t_a >= 0.0 && t_a <= 1.0 && boundary_closed?(m_val, lo, hi, inc_lo, inc_hi)
            builder.add_point(point_at_t(p1, p2, t_a, m1, m2, lo, hi))
          else
            builder.break!
          end
        end
      end

      # The point at parameter +t+ of the segment, with X/Y (and Z when
      # present) interpolated and M snapped inside [lo, hi]. Endpoints
      # are returned as-is so original vertices are preserved.
      def point_at_t(p1, p2, t, m1, m2, lo, hi)
        return p1 if t == 0.0
        return p2 if t == 1.0
        point = Point.new(srid, with_z, with_m)
        point.set_x_y(p1.x + t * (p2.x - p1.x), p1.y + t * (p2.y - p1.y))
        point.z = p1.z + t * (p2.z - p1.z) if with_z && p1.z && p2.z
        m_val = m1 + t * (m2 - m1)
        m_val = lo if m_val < lo
        m_val = hi if m_val > hi
        point.m = m_val
        point
      end

      def measure_covered?(v, lo, hi, inc_lo, inc_hi)
        (v > lo || (inc_lo && v == lo)) && (v < hi || (inc_hi && v == hi))
      end

      def boundary_closed?(v, lo, hi, inc_lo, inc_hi)
        (v == lo && inc_lo) || (v == hi && inc_hi)
      end

      # Accumulates isolated points and line pieces while segments are
      # scanned, merging pieces that meet on a junction point, and
      # builds the final geometry whose type reflects what was found.
      class MeasureResultBuilder
        def initialize(srid, with_z, with_m)
          @srid = srid
          @with_z = with_z
          @with_m = with_m
          @points = []
          @lines = []
          @run = nil
          @last = nil
        end

        # An isolated location. Duplicates produced twice at the junction
        # of two consecutive segments are reported only once.
        def add_point(point)
          return if @run && @run.last == point
          return if @run.nil? && @last == point
          flush_run
          @points << point
          @last = point
        end

        # A piece of line. Merged with the pending piece when they share
        # a junction point; a pending isolated point equal to the piece
        # start is absorbed by the piece.
        def add_piece(pa, pb)
          if @run
            if @run.last == pa
              @run << pb
              @last = pb
              return
            end
            flush_run
          end
          @points.pop if @points.last == pa
          @run = [pa, pb]
          @last = pb
        end

        # Cuts the continuity of the line (missing measure or segment
        # outside the interval): pending pieces and locations are closed.
        def break!
          flush_run
          @last = nil
        end

        def result
          flush_run
          return nil if @points.empty? && @lines.empty?
          return point_result if @lines.empty?
          return line_result if @points.empty?
          GeometryCollection.from_geometries(@points + @lines, @srid, @with_z, @with_m)
        end

        private

        def point_result
          return @points.first if @points.length == 1
          MultiPoint.from_points(@points, @srid, @with_z, @with_m)
        end

        def line_result
          return @lines.first if @lines.length == 1
          multi = MultiLineString.from_line_strings(@lines, @srid, @with_z, @with_m)
          multi.with_z = @with_z
          multi.with_m = @with_m
          multi
        end

        def flush_run
          return unless @run
          @lines << LineString.from_points(@run, @srid, @with_z, @with_m)
          @run = nil
        end
      end
    end

    class LineString
      include LinearReferencing
    end

    class MultiLineString
      include LinearReferencing
    end
  end
end
