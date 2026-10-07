# -*- coding: utf-8 -*-
require 'geo_ruby/simple_features/ellipsoid'

module GeoRuby
  module SimpleFeatures
    # Geodesic (shortest path) computations on an ellipsoid of
    # revolution, in the spirit of Karney's algorithms (GeographicLib)
    # but with the distance and longitude integrals evaluated through
    # Fourier series whose coefficients are computed numerically by the
    # trapezoid rule instead of truncated expansions in the flattening.
    #
    # The geodesic is mapped to a great circle on Bessel's auxiliary
    # sphere.  With sigma the arc length on that circle (measured from
    # the node), alpha0 the azimuth at the node and k = e' cos(alpha0),
    # the exact relations used are:
    #
    #   ds/dsigma = b * sqrt(1 + k^2 sin^2 sigma)
    #   dlambda   = domega - f sin(alpha0) (2 - f) /
    #               (1 + (1 - f) sqrt(1 + k^2 sin^2 sigma)) dsigma
    #
    # Both integrands are analytic and pi-periodic, so the trapezoid
    # rule used to build their Fourier series converges spectrally and
    # the series are accurate to machine precision.  The inverse
    # problem is solved by a safeguarded secant iteration on the
    # departure azimuth (with an inner fixed-point solve for the
    # spherical longitude difference), restarted from several
    # azimuths, so that near-antipodal configurations - where plain
    # Vincenty-style iteration fails to converge - are still solved
    # on the ellipsoid.  A spherical approximation is NEVER
    # substituted for the ellipsoidal result: when the iteration does
    # not converge, the returned result is flagged
    # <tt>converged? == false</tt> and carries an error estimate.
    module Geodesic
      DEG2RAD = Math::PI / 180.0
      RAD2DEG = 180.0 / Math::PI
      TWO_PI  = 2.0 * Math::PI

      # Fourier series of the distance and longitude integrands for
      # one geodesic (fixed equatorial azimuth alpha0).
      class Line
        HARMONICS = 12

        # Trapezoid points used to build the Fourier series.  The
        # trapezoid rule converges spectrally for these analytic
        # integrands, but the convergence rate degrades when k^2 gets
        # large (extreme flattening), so the resolution adapts to k^2.
        # WGS84 (k^2 <= 0.0068) always uses 64 points, which already
        # keeps the aliasing error far below machine precision.
        def self.quadrature_points(k2)
          if k2 > 100.0 then 4096
          elsif k2 > 10.0 then 1024
          elsif k2 > 1.0 then 256
          elsif k2 > 0.1 then 128
          else 64
          end
        end

        # Mean value of ds/dsigma / b (i.e. A1).
        attr_reader :distance_mean
        # Estimated truncation error of both series, in meters.
        attr_reader :truncation_error

        def initialize(ellipsoid, sin_alpha0, cos_alpha0)
          @b = ellipsoid.b
          @k2 = ellipsoid.ep2 * cos_alpha0 * cos_alpha0
          one_minus_f = 1.0 - ellipsoid.f
          k2 = @k2
          g = ->(sig) { Math.sqrt(1.0 + k2 * Math.sin(sig)**2) }
          h = lambda do |sig|
            (2.0 - ellipsoid.f) /
              (1.0 + one_minus_f * Math.sqrt(1.0 + k2 * Math.sin(sig)**2))
          end
          @dist = fourier_series(g)
          @lon  = fourier_series(h)
          @distance_mean = @dist[0]
          @truncation_error =
            @b * (tail_bound(@dist[1]) +
                  ellipsoid.f.abs * tail_bound(@lon[1]))
        end

        # Integral of ds/dsigma / b from 0 to sigma:
        #   A1 * sigma + sum C1[l] sin(2 l sigma)
        def distance_integral(sigma)
          series_eval(@dist, sigma)
        end

        # Exact derivative of distance_integral: sqrt(1 + k^2 sin^2).
        def distance_rate(sigma)
          Math.sqrt(1.0 + @k2 * Math.sin(sigma)**2)
        end

        # Integral of the longitude integrand from 0 to sigma.
        def longitude_integral(sigma)
          series_eval(@lon, sigma)
        end

        private

        # Decomposes a pi-periodic analytic function p into
        #   p(sigma) = a0 + sum_l 2 l c_l cos(2 l sigma)
        # and returns [a0, [c1, ..., cL]] so that the integral of p is
        #   P(sigma) = a0 * sigma + sum_l c_l sin(2 l sigma).
        def fourier_series(func)
          n = Line.quadrature_points(@k2)
          values = Array.new(n) { |j| func.call(Math::PI * j / n) }
          a0 = values.inject(0.0, :+) / n
          coeffs = (1..HARMONICS).map do |l|
            sum = 0.0
            n.times { |j| sum += values[j] * Math.cos(2.0 * l * Math::PI * j / n) }
            sum / (n * l)
          end
          [a0, coeffs]
        end

        def series_eval(series, sigma)
          a0, coeffs = series
          sum = a0 * sigma
          coeffs.each_with_index do |c, i|
            sum += c * Math.sin(2.0 * (i + 1) * sigma)
          end
          sum
        end

        # Rough bound on the dropped tail of the series, assuming
        # geometric decay of the coefficients.
       def tail_bound(coeffs)
         last = coeffs[-1].abs
         prev = coeffs[-2].abs
         ratio = prev.positive? ? [last / prev, 0.9].min : 0.0
         last / (1.0 - ratio)
       end
     end
      # Result of the inverse geodesic problem (see Point#geodesic_inverse).
      class InverseResult
        # Ellipsoidal distance between the two points, in meters.
        # Well defined even when the azimuths are not unique.
        attr_reader :distance
        # Forward azimuth at the first point, in degrees in [0, 360).
        # +nil+ when the two points coincide (azimuth meaningless).
        attr_reader :azimuth1
        # Forward azimuth at the second point, in degrees in [0, 360).
        # +nil+ when the two points coincide.
        attr_reader :azimuth2
        # Estimated upper bound of the error of +distance+, in meters.
        # It accounts for the truncation of the Fourier series and,
        # when the iteration did not fully converge, for the remaining
        # residual of the solver.
        attr_reader :error_estimate
        # Number of solver iterations (diagnostic).
        attr_reader :iterations

        def initialize(distance:, azimuth1:, azimuth2:, converged:,
                       error_estimate:, azimuths_unique:, iterations:)
          @distance = distance
          @azimuth1 = azimuth1
          @azimuth2 = azimuth2
          @converged = converged
          @error_estimate = error_estimate
          @azimuths_unique = azimuths_unique
          @iterations = iterations
        end

        # True when the solver converged to the requested tolerance.
        # When false, +distance+ is the best estimate found and
        # +error_estimate+ reflects the remaining residual; the result
        # is still computed on the ellipsoid (no spherical fallback).
        def converged?
          @converged
        end

        # True when the shortest geodesic - and therefore both
        # azimuths - is uniquely defined.  False for coincident points,
        # exactly antipodal configurations and pole-to-pole lines,
        # where several shortest geodesics exist; in those cases the
        # distance is still exact but the returned azimuths (when not
        # nil) are only one conventional choice.
        def azimuths_unique?
          @azimuths_unique
        end
      end

      # Result of the direct geodesic problem (see Point#geodesic_direct).
      class DirectResult
        # The destination Point (same SRID as the origin).
        attr_reader :point
        # Forward azimuth at the destination, in degrees in [0, 360).
        attr_reader :azimuth
        # Estimated upper bound of the position error, in meters.
        attr_reader :error_estimate
        # Number of solver iterations (diagnostic).
        attr_reader :iterations

        def initialize(point:, azimuth:, converged:, error_estimate:, iterations:)
          @point = point
          @azimuth = azimuth
          @converged = converged
          @error_estimate = error_estimate
          @iterations = iterations
        end

        def converged?
          @converged
        end
      end

      class << self
        # Inverse geodesic problem: distance and forward azimuths
        # between two Points (x = lon, y = lat, in degrees).
        # +ellipsoid+ may be an Ellipsoid, [a, f], {a:, f:} or nil
        # (defaults to WGS84).
        def inverse(point1, point2, ellipsoid = nil)
          ell = resolve_ellipsoid(ellipsoid)
          lat1 = check_latitude(point1.y)
          lon1 = normalize_longitude(point1.x)
          lat2 = check_latitude(point2.y)
          lon2 = normalize_longitude(point2.x)
          lon12 = lon2 - lon1
          lon12 -= 360.0 if lon12 > 180.0
          lon12 += 360.0 if lon12 <= -180.0
          inverse_degrees(ell, lat1, lon12, lat2)
        end

        # Direct geodesic problem: from +point+, travel +distance+
        # meters along the geodesic whose initial azimuth (clockwise
        # from north, in degrees) is +azimuth+.
        def direct(point, azimuth, distance, ellipsoid = nil)
          ell = resolve_ellipsoid(ellipsoid)
          lat1 = check_latitude(point.y)
          lon1 = normalize_longitude(point.x)
          azi1 = normalize_azimuth(azimuth)
          s12 = check_finite(distance, 'distance')
          if s12.negative?
            s12 = -s12
            azi1 = (azi1 + 180.0) % 360.0
          end
          lat2, lon2, azi2, err, its =
            direct_radians(ell, lat1 * DEG2RAD, lon1 * DEG2RAD,
                           azi1 * DEG2RAD, s12)
          DirectResult.new(
            point: Point.from_x_y(lon2 * RAD2DEG, lat2 * RAD2DEG, point.srid),
            azimuth: azi2 * RAD2DEG % 360.0,
            converged: true,
            error_estimate: err,
            iterations: its
          )
        end

        private

        # ----------------------------------------------------------
        # input handling
        # ----------------------------------------------------------

        def resolve_ellipsoid(ellipsoid)
          case ellipsoid
          when nil then Ellipsoid::WGS84
          when Ellipsoid then ellipsoid
          when Array then Ellipsoid.new(ellipsoid[0], ellipsoid[1])
          when Hash
            a = ellipsoid[:a] || ellipsoid['a']
            f = ellipsoid[:f] || ellipsoid['f']
            Ellipsoid.new(a, f)
          else
            raise ArgumentError,
                  "unsupported ellipsoid: #{ellipsoid.inspect}"
          end
        end

        def check_finite(value, name)
          unless value.is_a?(Numeric) && value.finite?
            raise ArgumentError,
                  "#{name} must be a finite number, got #{value.inspect}"
          end
          value.to_f
        end

        def check_latitude(lat)
          lat = check_finite(lat, 'latitude')
          if lat.abs > 90.0
            if lat.abs <= 90.0 + 1e-9
              return Math.copysign(90.0, lat)
            end
            raise ArgumentError,
                  "latitude out of range [-90, 90]: #{lat}"
          end
          lat
        end

        # Longitude in degrees, normalized to [-180, 180).
        def normalize_longitude(lon)
          lon = check_finite(lon, 'longitude') % 360.0
          lon -= 360.0 if lon >= 180.0
          lon
        end

        # Azimuth in degrees, normalized to [0, 360).
        def normalize_azimuth(azi)
          check_finite(azi, 'azimuth') % 360.0
        end

        # Longitude in radians, wrapped to [-pi, pi).
        def wrap_longitude(lam)
          lam %= TWO_PI
          lam -= TWO_PI if lam >= Math::PI
          lam
        end

        # ----------------------------------------------------------
        # direct problem
        # ----------------------------------------------------------

        # All angles in radians.  Returns [lat2, lon2, azi2, error, its].
        def direct_radians(ell, phi1, lam1, alpha1, s12)
          return [phi1, lam1, alpha1, 0.0, 0] if s12.zero?

          f = ell.f
          one_minus_f = 1.0 - f
          sb1 = Math.sin(phi1)
          cb1 = Math.cos(phi1)
          sa1 = Math.sin(alpha1)
          ca1 = Math.cos(alpha1)

          # Equatorial geodesic: exactly on the equator heading due
          # east or west (singular case of Bessel's formulation).
          # (cos/sin of pi/2 and pi are ~1e-17, not exactly zero.)
          if sb1.zero? && ca1.abs < 1e-15
            lam2 = lam1 + sa1 * s12 / ell.a
            return [0.0, wrap_longitude(lam2), alpha1, 0.0, 0]
          end

          # Meridional geodesic (due north/south), including
          # departures from a pole, where every geodesic is a meridian.
          if sa1.abs < 1e-15 || cb1.abs < 1e-15
            return direct_meridian(ell, phi1, lam1, alpha1, s12)
          end

          beta1 = Math.atan2(one_minus_f * sb1, cb1)
          sbet1 = Math.sin(beta1)
          cbet1 = Math.cos(beta1)
          sa0 = cbet1 * sa1
          ca0 = Math.hypot(ca1, sa1 * sbet1)
          sig1 = Math.atan2(sbet1, cbet1 * ca1)
          omg1 = Math.atan2(sa0 * Math.sin(sig1), Math.cos(sig1))
          line = Line.new(ell, sa0, ca0)

          # Invert the distance series by Newton's method (the series
          # is monotonic, so this converges from the first guess).
          target = s12 / ell.b
          base = line.distance_integral(sig1)
          sig2 = sig1 + target / line.distance_mean
          iterations = 0
          loop do
            delta = (line.distance_integral(sig2) - base - target) /
                    line.distance_rate(sig2)
            sig2 -= delta
            iterations += 1
            break if delta.abs < 1e-14 || iterations >= 50
          end

          ss2 = Math.sin(sig2)
          cs2 = Math.cos(sig2)
          sbet2 = ca0 * ss2
          cbet2 = Math.hypot(sa0, ca0 * cs2)
          lat2 = Math.atan2(sbet2, one_minus_f * cbet2)
          lam12 = (omega_unwrapped(sig2, sa0) - omg1) -
                  f * sa0 * (line.longitude_integral(sig2) -
                             line.longitude_integral(sig1))
          lon2 = wrap_longitude(lam1 + lam12)
          azi2 = Math.atan2(sa0, ca0 * cs2)
          [lat2, lon2, azi2, line.truncation_error, iterations]
        end

        # Meridional direct problem.  Handles due north/south
        # departures and departures from an exact pole.  Convention at
        # the poles (azimuth is meaningless there without one): the
        # geodesic leaves the north pole along the meridian
        # lam1 + pi - alpha1 and the south pole along lam1 + alpha1.
        def direct_meridian(ell, phi1, lam1, alpha1, s12)
          f = ell.f
          line = Line.new(ell, 0.0, 1.0) # k = e' (meridian)
          if phi1.abs > Math::PI / 2 - 1e-12
            north_pole = phi1.positive?
            lam_d = north_pole ? lam1 + Math::PI - alpha1 : lam1 + alpha1
            sig1 = phi1
            travel_north = !north_pole
          else
            lam_d = lam1
            sig1 = Math.atan2((1.0 - f) * Math.sin(phi1), Math.cos(phi1))
            # in the meridian branch alpha1 is ~0 (north) or ~pi (south)
            travel_north = alpha1 < Math::PI / 2
          end
          goal = line.distance_integral(sig1) +
                 (travel_north ? 1.0 : -1.0) * s12 / ell.b
          sig2 = sig1 + (travel_north ? 1.0 : -1.0) * s12 / ell.b /
                line.distance_mean
          iterations = 0
          loop do
            delta = (line.distance_integral(sig2) - goal) /
                    line.distance_rate(sig2)
            sig2 -= delta
            iterations += 1
            break if delta.abs < 1e-14 || iterations >= 50
          end
          # Wrap sigma into [-pi/2, pi/2], tracking pole crossings.
          # The tiny bias makes an exact landing on a pole count as
          # 'arriving' rather than 'crossing'.
          crossings = ((sig2 + Math::PI / 2 - 1e-13) / Math::PI).floor
          beta2 =
            if crossings.odd?
              crossings * Math::PI - sig2
            else
              sig2 - crossings * Math::PI
            end
          base_azi = travel_north ? 0.0 : Math::PI
          azi2 = crossings.odd? ? Math::PI - base_azi : base_azi
          lat2 = Math.atan2(Math.sin(beta2), (1.0 - f) * Math.cos(beta2))
          lon2 = wrap_longitude(lam_d + crossings * Math::PI)
          [lat2, lon2, azi2, line.truncation_error, iterations]
        end

        # Spherical longitude on the auxiliary circle, unwrapped so
        # that it grows continuously with sigma.
       def omega_unwrapped(sig, sa0)
         principal = Math.atan2(Math.sin(sig), Math.cos(sig))
         Math.atan2(sa0 * Math.sin(sig), Math.cos(sig)) +
           TWO_PI * ((sig - principal) / TWO_PI).round
       end

        # ----------------------------------------------------------
        # inverse problem
        # ----------------------------------------------------------

        # lat/lon differences in degrees; lon12 in (-180, 180].
        def inverse_degrees(ell, lat1d, lon12d, lat2d)
          # Coincident points (two poles with the same sign are the
          # same physical point): the distance is exactly zero and
          # the azimuths are undefined.
          if (lat1d == lat2d && lon12d.zero?) ||
             (lat1d.abs == 90.0 && lat1d == lat2d)
            return InverseResult.new(
              distance: 0.0, azimuth1: nil, azimuth2: nil,
              converged: true, error_estimate: 0.0,
              azimuths_unique: false, iterations: 0
            )
          end

          phi1 = lat1d * DEG2RAD
          phi2 = lat2d * DEG2RAD
          lam12 = lon12d * DEG2RAD

          if lat1d.abs == 90.0 || lat2d.abs == 90.0
            return inverse_with_pole(ell, phi1, lam12, phi2)
          end

          beta1 = reduced_latitude(ell, phi1)
          beta2 = reduced_latitude(ell, phi2)

          # Meridional cases.  Longitude differences within 1e-9 rad of
          # 0 or pi are treated as exactly meridional: the true geodesic
          # is then indistinguishable from the meridian one (the
          # distance differs by far less than 1e-9 m and the azimuths
          # by less than 1e-5 deg).
          if lam12.abs <= 1e-9
            # Same meridian.
            line = Line.new(ell, 0.0, 1.0)
            s12 = ell.b * (line.distance_integral(beta2) -
                           line.distance_integral(beta1)).abs
            azi = beta2 > beta1 ? 0.0 : 180.0
            return InverseResult.new(
              distance: s12, azimuth1: azi, azimuth2: azi,
              converged: true,
              error_estimate: line.truncation_error,
              azimuths_unique: true, iterations: 0
            )
          end

          if (lam12.abs - Math::PI).abs <= 1e-9
            return inverse_opposite_meridians(ell, beta1, beta2)
          end

          # Mirror the longitudes so that lam12 > 0; the azimuths are
          # un-mirrored at the end.
          mirror = lam12.negative?
          lam = lam12.abs

          if phi1.zero? && phi2.zero?
            result = inverse_equatorial(ell, lam)
          else
            result = inverse_general(ell, phi1, phi2, lam)
          end

          azi1 = result[:azi1] * RAD2DEG
          azi2 = result[:azi2] * RAD2DEG
          if mirror
            azi1 = 360.0 - azi1
            azi2 = 360.0 - azi2
          end
          InverseResult.new(
            distance: result[:s12],
            azimuth1: azi1 % 360.0,
            azimuth2: azi2 % 360.0,
            converged: result[:converged],
            error_estimate: result[:error],
            azimuths_unique: true,
            iterations: result[:iterations]
          )
        end

        def reduced_latitude(ell, phi)
          Math.atan2((1.0 - ell.f) * Math.sin(phi), Math.cos(phi))
        end

        # Inverse problem when at least one point is an exact pole.
        # Every geodesic through a pole is a meridian, so the distance
        # is exact; the azimuth at the pole is only a convention (see
        # direct_meridian), hence azimuths_unique? is false.
        def inverse_with_pole(ell, phi1, lam12, phi2)
          line = Line.new(ell, 0.0, 1.0)
          quadrant = line.distance_integral(Math::PI / 2)
          err = line.truncation_error
          pole1 = phi1.abs > Math::PI / 2 - 1e-12
          pole2 = phi2.abs > Math::PI / 2 - 1e-12
          if pole1 && pole2
            # Pole to pole: infinitely many shortest meridians.
            return InverseResult.new(
              distance: 2.0 * ell.b * quadrant,
              azimuth1: nil, azimuth2: nil, converged: true,
              error_estimate: err, azimuths_unique: false, iterations: 0
            )
          end
          if pole1
            beta2 = reduced_latitude(ell, phi2)
            arc = quadrant - (phi1.positive? ? 1.0 : -1.0) *
                  line.distance_integral(beta2)
            # Departure azimuth follows the pole convention; the
            # arrival azimuth is due south (from N pole) or due north
            # (from S pole).
            azi1 = phi1.positive? ? Math::PI - lam12 : lam12
            azi2 = phi1.positive? ? Math::PI : 0.0
          else
            beta1 = reduced_latitude(ell, phi1)
            # arc from point 1 to the pole at phi2:
            #   N pole: S(pi/2) - S(beta1);  S pole: S(beta1) - S(-pi/2)
            arc = quadrant - (phi2.positive? ? 1.0 : -1.0) *
                  line.distance_integral(beta1)
            azi1 = phi2.positive? ? 0.0 : Math::PI
            azi2 = phi2.positive? ? 0.0 : Math::PI
          end
          InverseResult.new(
            distance: ell.b * arc,
            azimuth1: (azi1 * RAD2DEG) % 360.0,
            azimuth2: (azi2 * RAD2DEG) % 360.0,
            converged: true, error_estimate: err,
            azimuths_unique: false, iterations: 0
          )
        end

        # Points on opposite meridians (lon12 = 180): the shortest
        # geodesic runs over one of the poles.  When beta1 + beta2 = 0
        # both polar paths have the same length and the azimuth is not
        # unique (the distance still is).
        def inverse_opposite_meridians(ell, beta1, beta2)
          line = Line.new(ell, 0.0, 1.0)
          quadrant = line.distance_integral(Math::PI / 2)
          arc_north = 2.0 * quadrant -
                      line.distance_integral(beta1) -
                      line.distance_integral(beta2)
          arc_south = 2.0 * quadrant +
                      line.distance_integral(beta1) +
                      line.distance_integral(beta2)
          if arc_north <= arc_south
            s12 = ell.b * arc_north
            azi1 = 0.0
            azi2 = 180.0
          else
            s12 = ell.b * arc_south
            azi1 = 180.0
            azi2 = 0.0
          end
          InverseResult.new(
            distance: s12, azimuth1: azi1, azimuth2: azi2,
            converged: true, error_estimate: line.truncation_error,
            azimuths_unique: arc_north != arc_south, iterations: 0
          )
        end

        # Both points exactly on the equator.  The equator itself is a
        # geodesic, but beyond the cut-off lon12 > pi (1 - f) the
        # shortest path leaves the equator; both candidates are
        # evaluated and the shorter one wins.
        def inverse_equatorial(ell, lam)
          best = { s12: ell.a * lam, azi1: Math::PI / 2,
                   azi2: Math::PI / 2, converged: true,
                   error: 0.0, iterations: 0 }
          if lam > Math::PI * (1.0 - ell.f) - 1e-9
            gen = inverse_general(ell, 0.0, 0.0, lam)
           best = gen if gen[:converged] && gen[:s12] < best[:s12]
         end
         best
       end

        # General inverse problem for lam in (0, pi), away from the
        # meridional and polar special cases.  Solves
        #   beta(sig2(alpha1)) = beta2
        # for the departure azimuth alpha1 with a secant iteration,
        # where sig2 follows from the spherical longitude omega2 and
        # omega12 is slaved to the longitude equation by a
        # fixed-point iteration.  Several starting azimuths (a
        # spherical guess plus a coverage set and a sign-change scan)
        # are used so that near-antipodal cases still converge on the
        # ellipsoid; the shortest converged geodesic is returned.
        def inverse_general(ell, phi1, phi2, lam)
          beta1 = reduced_latitude(ell, phi1)
          beta2 = reduced_latitude(ell, phi2)
          sbet1 = Math.sin(beta1)
          cbet1 = Math.cos(beta1)
          sbet2 = Math.sin(beta2)
          cbet2 = Math.cos(beta2)

          eval_g = build_residual(ell, sbet1, cbet1, beta2, lam)

          guess = Math.atan2(Math.sin(lam) * cbet2,
                             cbet1 * sbet2 - sbet1 * cbet2 * Math.cos(lam))
          # Fast path: away from the antipodal region the spherical
          # guess lands in the basin of the shortest geodesic.
          if lam <= 2.5
            cand = secant_inverse(ell, eval_g, guess)
            return cand if cand[:converged]
          end

          # Sign-change brackets from a coarse scan; each bracket is
          # refined by bisection, which cannot stall in the flat
          # near-antipodal region of the residual the way the secant
          # method can.
          brackets = []
          prev_x = 1e-9
          prev_g = eval_g.call(prev_x)
          prev_g = prev_g && prev_g[:g]
          24.times do |i|
            x = (i + 1).to_f / 24 * (Math::PI - 2e-9) + 1e-9
            r = eval_g.call(x)
            g = r && r[:g]
            brackets << [prev_x, x] if g && prev_g && g * prev_g < 0
            if g
              prev_x = x
              prev_g = g
            end
          end

          candidates = []
          iterations = 0
          brackets.each do |a, b|
            cand = bisect_inverse(ell, eval_g, a, b)
            iterations += cand[:iterations]
            candidates << cand if cand[:converged]
          end

          # Secant from a coverage set of starting azimuths, in case a
          # root has no sign change (tangent residual).  When brackets
          # were found, the bisections already cover the crossing
          # roots and two starts suffice.
          starts = [guess, Math::PI / 2]
          if brackets.empty?
            starts += [Math::PI / 6, 5.0 * Math::PI / 6,
                       Math::PI / 12, 11.0 * Math::PI / 12,
                       Math::PI / 4, 3.0 * Math::PI / 4]
          end
          starts.each do |start|
            cand = secant_inverse(ell, eval_g, start)
            iterations += cand[:iterations]
            candidates << cand if cand[:converged]
          end

          # The shortest geodesic always has sig12 <= pi; longer
          # (still valid) geodesics are only a fallback.
          short = candidates.select { |c| c[:sig12] <= Math::PI + 1e-9 }
          pool = short.empty? ? candidates : short
          best = pool.min_by { |c| c[:s12] }
          best ||= secant_inverse(ell, eval_g, guess)
          best[:iterations] = iterations
          best
        end

        # Bisection on a sign-change bracket: slow but unconditionally
        # convergent, used for the hard (near-antipodal) cases.
        def bisect_inverse(ell, eval_g, a, b)
          ra = eval_g.call(a)
          rb = eval_g.call(b)
          iterations = 2
          if ra.nil? || rb.nil?
            return { converged: false, s12: 0.0, azi1: a, azi2: 0.0,
                     sig12: 0.0, residual: Float::INFINITY,
                     iterations: iterations, error: Float::INFINITY }
          end
          return finalize_candidate(ell, ra, a, true, iterations) if ra[:g].zero?
          return finalize_candidate(ell, rb, b, true, iterations) if rb[:g].zero?
          ga = ra[:g]
          60.times do
            m = (a + b) / 2
            rm = eval_g.call(m)
            iterations += 1
            next if rm.nil?
            if rm[:g].abs < 1e-13
              return finalize_candidate(ell, rm, m, true, iterations)
            end
            if ga * rm[:g] < 0
              b = m
            else
              a = m
              ga = rm[:g]
            end
          end
          m = (a + b) / 2
          rm = eval_g.call(m)
          if rm.nil?
            return { converged: false, s12: 0.0, azi1: m, azi2: 0.0,
                     sig12: 0.0, residual: Float::INFINITY,
                     iterations: iterations + 1, error: Float::INFINITY }
          end
          finalize_candidate(ell, rm, m, rm[:g].abs < 1e-13, iterations + 1)
        end

        # Builds the residual function G(alpha1) = beta2_computed -
        # beta2 (radians).  Each evaluation also returns the resolved
        # geometry of the trial geodesic.
        def build_residual(ell, sbet1, cbet1, beta2, lam)
          f = ell.f
          lambda do |x|
            sa1 = Math.sin(x)
            ca1 = Math.cos(x)
            sa0 = cbet1 * sa1
            return nil if sa0.zero?
            ca0 = Math.hypot(ca1, sa1 * sbet1)
            sig1 = Math.atan2(sbet1, cbet1 * ca1)
            omg1 = Math.atan2(sa0 * Math.sin(sig1), Math.cos(sig1))
            line = Line.new(ell, sa0, ca0)
            omg12 = lam
            sig12 = 0.0
            50.times do
              omg2 = omg1 + omg12
              sig2p = Math.atan2(Math.sin(omg2), sa0 * Math.cos(omg2))
              sig12 = sig2p - sig1
              sig12 += TWO_PI if sig12.negative?
              delta = lam + f * sa0 *
                      (line.longitude_integral(sig1 + sig12) -
                       line.longitude_integral(sig1)) - omg12
              omg12 += delta
              break if delta.abs < 1e-16
            end
            sbet2c = ca0 * Math.sin(sig1 + sig12)
            sbet2c = 1.0 if sbet2c > 1.0
            sbet2c = -1.0 if sbet2c < -1.0
            { g: Math.asin(sbet2c) - beta2, sig1: sig1, sig12: sig12,
              line: line, sa0: sa0, ca0: ca0 }
          end
        end

        # Secant iteration on the departure azimuth from one starting
        # value.  Always ellipsoidal: on failure it reports
        # converged: false together with the best estimate found.
        def secant_inverse(ell, eval_g, start)
          x0 = start.clamp(1e-9, Math::PI - 1e-9)
          x1 = x0 + (x0 < Math::PI / 2 ? 0.01 : -0.01)
          r0 = eval_g.call(x0)
          r1 = eval_g.call(x1)
          iterations = 2
          best = nil
          best_x = nil
          [[x0, r0], [x1, r1]].each do |x, r|
            next unless r
            if best.nil? || r[:g].abs < best[:g].abs
              best = r
              best_x = x
            end
          end
          20.times do
            break if r0.nil? || r1.nil?
            denom = r1[:g] - r0[:g]
            break if denom.zero?
            x2 = x1 - r1[:g] * (x1 - x0) / denom
            break if x2.nan? || x2.infinite?
            x2 = x2.clamp(1e-9, Math::PI - 1e-9)
            r2 = eval_g.call(x2)
            iterations += 1
            break if r2.nil?
            if best.nil? || r2[:g].abs < best[:g].abs
              best = r2
              best_x = x2
            end
            break if r2[:g].abs < 1e-13
            x0 = x1
            r0 = r1
            x1 = x2
            r1 = r2
          end
          if best.nil?
            return { converged: false, s12: 0.0, azi1: start, azi2: 0.0,
                     sig12: 0.0, residual: Float::INFINITY,
                     iterations: iterations, error: Float::INFINITY }
          end
          converged = best[:g].abs < 1e-13
          finalize_candidate(ell, best, best_x, converged, iterations)
        end

        def finalize_candidate(ell, r, x, converged, iterations)
          sig1 = r[:sig1]
          sig2 = sig1 + r[:sig12]
          line = r[:line]
          s12 = ell.b * (line.distance_integral(sig2) -
                         line.distance_integral(sig1))
          azi2 = Math.atan2(r[:sa0], r[:ca0] * Math.cos(sig2))
          { converged: converged, s12: s12, azi1: x, azi2: azi2,
            sig12: r[:sig12],
            residual: r[:g].abs, iterations: iterations,
            error: line.truncation_error + ell.b * r[:g].abs }
        end
      end
    end
  end
end
