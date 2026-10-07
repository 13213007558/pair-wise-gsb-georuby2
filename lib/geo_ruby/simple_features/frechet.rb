require 'geo_ruby/simple_features/geometry'

module GeoRuby
  module SimpleFeatures
    # Continuous Frechet distance between polygonal curves in the plane,
    # after Alt & Godau (1995).
    #
    # A matching is a pair of monotone reparameterizations of the two
    # curves: both parameters may only advance (or stay put) along their
    # own curve. The continuous Frechet distance is the smallest epsilon
    # such that a matching exists keeping paired points within epsilon.
    # This is computed exactly on the free-space diagram (one cell per
    # pair of segments) in O(m*n) per decision, not by enumerating
    # discrete vertex alignments.
    module Frechet
      # Tolerance used when testing whether an interval endpoint touches
      # a cell corner (coordinates 0 and 1 in cell-local space).
      CORNER_TOL = 1e-12

      # Result of a distance computation. The true continuous Frechet
      # distance is guaranteed to lie in [lower_bound, upper_bound].
      class Result
        attr_reader :lower_bound, :upper_bound

        def initialize(lower_bound, upper_bound)
          @lower_bound = lower_bound
          @upper_bound = upper_bound
        end

        # Midpoint of the certified interval.
        def distance
          0.5 * (@lower_bound + @upper_bound)
        end

        # Width of the certified interval.
        def error_bound
          @upper_bound - @lower_bound
        end

        def to_s
          "Frechet distance in [#{@lower_bound}, #{@upper_bound}]"
        end
      end

      # Outcome of the decision problem "is the continuous Frechet
      # distance <= epsilon?". When the answer is positive, +witness+
      # holds a monotone matching as a list of [s, t] parameter pairs
      # (s along the first curve, t along the second, both in curve
      # parameter space where integer values are vertices). The witness
      # can be re-checked independently with Frechet.valid_witness?.
      class Decision
        attr_reader :epsilon, :witness

        def initialize(within, epsilon, witness = nil)
          @within = within
          @epsilon = epsilon
          @witness = witness
        end

        def within?
          @within
        end

        def to_s
          "Frechet distance #{within? ? '<=' : '>'} #{@epsilon}"
        end
      end

      class << self
        # Decision problem: is delta_F(pa, qa) <= epsilon?
        # Curves are arrays of [x, y] pairs. Free space is closed, so a
        # threshold that exactly touches the free-space boundary counts
        # as reachable.
        def decision(pa, qa, epsilon)
          check_curves(pa, qa)
          epsilon = epsilon.to_f
          eps2 = epsilon * epsilon

          if pa.length == 1 || qa.length == 1
            within, witness = degenerate_decision(pa, qa, eps2)
            return Decision.new(within, epsilon, witness)
          end

          return Decision.new(false, epsilon) unless endpoints_free?(pa, qa, eps2)

          fv, fh = free_edges(pa, qa, eps2)
          rv, rh = reachable_edges(pa.length, qa.length, fv, fh)
          within = !(rv[-1][-1].nil? && rh[-1][-1].nil?)
          witness = within ? build_witness(pa.length, qa.length, rv, rh) : nil
          Decision.new(within, epsilon, witness)
        end

        # Continuous Frechet distance with a certified error bound:
        # the exact value lies in [result.lower_bound, result.upper_bound],
        # and upper_bound - lower_bound <= tolerance * max(1, distance).
        def distance(pa, qa, tolerance = 1e-9)
          check_curves(pa, qa)

          if pa.length == 1 || qa.length == 1
            d = degenerate_distance(pa, qa)
            return Result.new(d, d)
          end

          lo = Math.sqrt([dist2(pa.first, qa.first), dist2(pa[-1], qa[-1])].max)
          hi = Math.sqrt([pa.map { |p| dist2(p, qa.first) }.max,
                          qa.map { |q| dist2(pa[-1], q) }.max].max)
          hi = lo if hi < lo
          tol_abs = tolerance * [1.0, lo].max

          200.times do
            break if hi - lo <= tol_abs
            mid = 0.5 * (lo + hi)
            if decision(pa, qa, mid).within?
              hi = mid
            else
              lo = mid
            end
          end
          Result.new(lo, hi)
        end

        # Independent checker for a monotone matching witness: verifies
        # endpoints, monotonicity, and that every matched pair (sampled
        # at all cell-boundary crossings, which suffices by convexity of
        # the free space inside each cell) lies within epsilon.
        def valid_witness?(pa, qa, witness, epsilon, tol = 1e-9)
          return false if witness.nil? || witness.empty?
          return false unless witness.first == [0.0, 0.0]
          return false unless witness.last == [(pa.length - 1).to_f, (qa.length - 1).to_f]

          limit = (epsilon + tol) * (epsilon + tol)
          witness.each_cons(2) do |a, b|
            return false if b[0] < a[0] - CORNER_TOL || b[1] < a[1] - CORNER_TOL

            crossing_params(a, b).each do |u|
              s = a[0] + u * (b[0] - a[0])
              t = a[1] + u * (b[1] - a[1])
              return false if dist2(point_at(pa, s), point_at(qa, t)) > limit
            end
          end
          true
        end

        # Point at curve parameter s (integer values are vertices).
        def point_at(curve, s)
          return curve[0] if curve.length == 1

          s = 0.0 if s.negative?
          s = (curve.length - 1).to_f if s > curve.length - 1
          k = s.floor
          k = curve.length - 2 if k >= curve.length - 1
          f = s - k
          ax, ay = curve[k]
          bx, by = curve[k + 1]
          [ax + f * (bx - ax), ay + f * (by - ay)]
        end

        private

        def check_curves(pa, qa)
          raise ArgumentError, 'empty curve' if pa.empty? || qa.empty?
        end

        def dist2(p, q)
          (p[0] - q[0])**2 + (p[1] - q[1])**2
        end

        def endpoints_free?(pa, qa, eps2)
          dist2(pa.first, qa.first) <= eps2 && dist2(pa[-1], qa[-1]) <= eps2
        end

        # Interval of t in [0, 1] with |a + t*(b-a) - p|^2 <= eps2,
        # or nil when empty. Handles zero-length segments (a == b).
        def free_interval(p, a, b, eps2)
          dx = b[0] - a[0]
          dy = b[1] - a[1]
          ex = a[0] - p[0]
          ey = a[1] - p[1]
          qa = dx * dx + dy * dy
          if qa.zero?
            return ex * ex + ey * ey <= eps2 ? [0.0, 1.0] : nil
          end

          qb = 2.0 * (dx * ex + dy * ey)
          qc = ex * ex + ey * ey - eps2
          disc = qb * qb - 4.0 * qa * qc
          return nil if disc.negative?

          root = Math.sqrt(disc)
          lo = (-qb - root) / (2.0 * qa)
          hi = (-qb + root) / (2.0 * qa)
          lo = 0.0 if lo.negative?
          hi = 1.0 if hi > 1.0
          lo > hi ? nil : [lo, hi]
        end

        # Free intervals on all cell edges of the diagram.
        # fv[i][j]: vertical edge, first curve fixed at vertex i,
        #           second curve moving along segment j.
        # fh[i][j]: horizontal edge, second curve fixed at vertex j,
        #           first curve moving along segment i.
        def free_edges(pa, qa, eps2)
          fv = Array.new(pa.length) do |i|
            Array.new(qa.length - 1) { |j| free_interval(pa[i], qa[j], qa[j + 1], eps2) }
          end
          fh = Array.new(pa.length - 1) do |i|
            Array.new(qa.length) { |j| free_interval(qa[j], pa[i], pa[i + 1], eps2) }
          end
          [fv, fh]
        end

        # interval intersected with [lo_bound, 1]
        def intersect_from(interval, lo_bound)
          return nil if interval.nil?

          lo = [interval[0], lo_bound].max
          lo <= interval[1] ? [lo, interval[1]] : nil
        end

        # Alt-Godau propagation of reachable intervals. rv[i][j] is the
        # reachable interval on the left edge of cell (i, j), rh[i][j]
        # on its bottom edge; each is nil or [lo, hi] in cell-local
        # coordinates.
        def reachable_edges(m, n, fv, fh)
          rv = Array.new(m - 1) { Array.new(n - 1) }
          rh = Array.new(m - 1) { Array.new(n - 1) }
          rv[0][0] = fv[0][0]
          rh[0][0] = fh[0][0]

          (0...(m - 1)).each do |i|
            (0...(n - 1)).each do |j|
              # propagation along the outer boundary of the diagram
              if i.zero? && j.positive? && rv[0][j].nil?
                prev = rv[0][j - 1]
                rv[0][j] = fv[0][j] if prev && prev[1] >= 1.0 - CORNER_TOL
              end
              if j.zero? && i.positive? && rh[i][0].nil?
                prev = rh[i - 1][0]
                rh[i][0] = fh[i][0] if prev && prev[1] >= 1.0 - CORNER_TOL
              end

              # a reachable corner makes both adjacent edges fully
              # reachable (free space inside a cell is convex)
              if (rv[i][j] && rv[i][j][0] <= CORNER_TOL) ||
                 (rh[i][j] && rh[i][j][0] <= CORNER_TOL)
                rv[i][j] = fv[i][j]
                rh[i][j] = fh[i][j]
              end

              if i + 1 < m - 1
                rv[i + 1][j] =
                  if rh[i][j]
                    fv[i + 1][j]
                  elsif rv[i][j]
                    intersect_from(fv[i + 1][j], rv[i][j][0])
                  end
              end
              if j + 1 < n - 1
                rh[i][j + 1] =
                  if rv[i][j]
                    fh[i][j + 1]
                  elsif rh[i][j]
                    intersect_from(fh[i][j + 1], rh[i][j][0])
                  end
              end
            end
          end
          [rv, rh]
        end

        # Backtrack a monotone path through the reachable cells, from
        # (m-1, n-1) back to (0, 0). Returned in forward order.
        def build_witness(m, n, rv, rh)
          i = m - 2
          j = n - 2
          s = (m - 1).to_f
          t = (n - 1).to_f
          path = [[s, t]]

          until s.zero? && t.zero?
            if s.zero?
              # on the left boundary the whole segment down to (0, 0)
              # is reachable; record the corners it crosses
              t = step_down(t)
              path << [0.0, t]
              next
            elsif t.zero?
              s = step_down(s)
              path << [s, 0.0]
              next
            end

            rvl = rv[i][j]
            rhl = rh[i][j]
            on_right = s == (i + 1).to_f
            on_top = t == (j + 1).to_f
            go_left =
              if on_right && on_top
                !rvl.nil?
              elsif on_right
                !rvl.nil? && t - j >= rvl[0] - CORNER_TOL
              else
                rhl.nil? || s - i < rhl[0] - CORNER_TOL
              end

            if go_left
              raise 'unreachable state in witness backtracking' if rvl.nil?

              tl = [[t - j, rvl[0]].max, rvl[1]].min
              t = j + tl
              s = i.to_f
              path << [s, t]
              i -= 1
            else
              raise 'unreachable state in witness backtracking' if rhl.nil?

              sl = [[s - i, rhl[0]].max, rhl[1]].min
              s = i + sl
              t = j.to_f
              path << [s, t]
              j -= 1
            end
          end
          path.reverse
        end

        # Next integer parameter value strictly below x (or 0).
        def step_down(x)
          f = x.floor.to_f
          nxt = f >= x ? x - 1.0 : f
          nxt.negative? ? 0.0 : nxt
        end

        # Parameters u in [0, 1] where the segment a->b in parameter
        # space crosses a cell boundary, plus the endpoints.
        def crossing_params(a, b)
          us = [0.0, 1.0]
          2.times do |axis|
            span = b[axis] - a[axis]
            next if span.zero?

            lo = [a[axis], b[axis]].min
            hi = [a[axis], b[axis]].max
            (lo.floor + 1).upto(hi.ceil) do |k|
              u = (k - a[axis]) / span
              us << u if u > 0.0 && u < 1.0
            end
          end
          us.uniq.sort
        end

        # Single-point curves: the matching pairs the fixed point with
        # every point of the other curve; the worst pair sits at a
        # vertex (distance to a fixed point is convex along a segment).
        def degenerate_decision(pa, qa, eps2)
          fixed, moving, flip = pa.length == 1 ? [pa[0], qa, false] : [qa[0], pa, true]
          within = moving.all? { |q| dist2(fixed, q) <= eps2 }
          witness = nil
          if within
            witness = (0...moving.length).map do |k|
              flip ? [k.to_f, 0.0] : [0.0, k.to_f]
            end
          end
          [within, witness]
        end

        def degenerate_distance(pa, qa)
          fixed, moving = pa.length == 1 ? [pa[0], qa] : [qa[0], pa]
          Math.sqrt(moving.map { |q| dist2(fixed, q) }.max)
        end
      end
    end

    class LineString < Geometry
      # Continuous Frechet distance to another line string, in the
      # plane of the x/y coordinates. Returns a Frechet::Result whose
      # [lower_bound, upper_bound] brackets the exact value.
      #
      # Closed line strings are compared as given, starting and ending
      # at points.first (no rotation to a "best" start point).
      def frechet_distance(other, tolerance = 1e-9)
        Frechet.distance(frechet_xy, other.send(:frechet_xy), tolerance)
      end

      # Decision problem: is the continuous Frechet distance to +other+
      # at most +epsilon+? Returns a Frechet::Decision; when positive it
      # carries a monotone matching witness for review.
      def frechet_decision(other, epsilon)
        Frechet.decision(frechet_xy, other.send(:frechet_xy), epsilon)
      end

      # Boolean shortcut for frechet_decision(other, epsilon).within?
      def frechet_within?(other, epsilon)
        frechet_decision(other, epsilon).within?
      end

      private

      def frechet_xy
        raise ArgumentError, 'cannot compute Frechet distance of an empty line string' if @points.empty?

        @points.map { |p| [p.x.to_f, p.y.to_f] }
      end
    end
  end
end
