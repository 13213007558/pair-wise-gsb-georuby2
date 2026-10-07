# frozen_string_literal: true

module GeoRuby
  module SimpleFeatures
    class Ellipsoid
      DEG2RAD = Math::PI / 180.0
      RAD2DEG = 180.0 / Math::PI

      attr_reader :semi_major_axis, :semi_minor_axis, :flattening

      def initialize(semi_major_axis:, semi_minor_axis:)
        @semi_major_axis = finite_positive_number!(semi_major_axis, :semi_major_axis)
        @semi_minor_axis = finite_positive_number!(semi_minor_axis, :semi_minor_axis)
        if @semi_minor_axis > @semi_major_axis
          raise ArgumentError, 'semi_minor_axis cannot be greater than semi_major_axis'
        end

        @flattening = (@semi_major_axis - @semi_minor_axis) / @semi_major_axis
      end

      def finite_positive_number!(value, name)
        unless value.is_a?(Numeric) && value.finite? && value > 0
          raise ArgumentError, "#{name} must be a finite positive number"
        end

        value.to_f
      end

      WGS84 = new(
        semi_major_axis: 6_378_137.0,
        semi_minor_axis: 6_356_752.314245179
      )
    end

    GeodesicSolution = Struct.new(
      :status,
      :converged,
      :distance,
      :initial_forward_azimuth,
      :final_forward_azimuth,
      :reverse_azimuth,
      :endpoint,
      :iterations,
      :distance_error,
      :longitude_error,
      :azimuth_error,
      :message,
      keyword_init: true
    ) do
      def converged?
        converged
      end

      def ambiguous_azimuth?
        status == :ambiguous_azimuth
      end

      def same_point?
        status == :same_point
      end
    end

    # Ellipsoidal geodesic calculations using Vincenty's auxiliary sphere.
    # The inverse problem uses fixed-point iteration for ordinary pairs and a
    # bracketed root when that iteration cannot make progress. It therefore does
    # not substitute a spherical answer near antipodal points.
    module Geodesic
      TOLERANCE = 1e-12
      MAX_ITERATIONS = 100

      module_function

      def inverse(point, other, ellipsoid = Ellipsoid::WGS84)
        ellipsoid = validate_ellipsoid(ellipsoid)
        point = validate_point!(point)
        other = validate_point!(other)
        lat1 = point.y * Ellipsoid::DEG2RAD
        lat2 = other.y * Ellipsoid::DEG2RAD
        lon_diff = wrap_degrees(other.x - point.x) * Ellipsoid::DEG2RAD
        sin_u1, cos_u1 = reduced_latitude(lat1, ellipsoid.flattening)
        sin_u2, cos_u2 = reduced_latitude(lat2, ellipsoid.flattening)

        at_pole = point.y.abs == 90.0 && other.y == point.y
        same_latitude = point.y == other.y
        same_position = at_pole || (same_latitude && lon_diff.zero?)
        if same_position
          return GeodesicSolution.new(
            status: :same_point,
            converged: true,
            distance: 0.0,
            initial_forward_azimuth: nil,
            final_forward_azimuth: nil,
            reverse_azimuth: nil,
            iterations: 0,
            distance_error: 0.0,
            longitude_error: 0.0,
            azimuth_error: 0.0,
            message: 'The points are identical; azimuth is undefined.'
          )
        end

        antipodal = point.y == -other.y && lon_diff.abs == Math::PI &&
                    point.y.abs != 90.0
        if antipodal
          distance = meridian_distance(-Math::PI / 2.0, Math::PI / 2.0, ellipsoid)
          return GeodesicSolution.new(
            status: :ambiguous_azimuth,
            converged: true,
            distance: distance,
            initial_forward_azimuth: nil,
            final_forward_azimuth: nil,
            reverse_azimuth: nil,
            iterations: 0,
            distance_error: 0.0,
            longitude_error: 0.0,
            azimuth_error: nil,
            message: 'The endpoints are antipodal; more than one shortest geodesic exists.'
          )
        end

        if lon_diff.zero?
          distance = meridian_distance(
            lat1,
          lat2,
          ellipsoid
          )
          initial = lat2 > lat1 ? 0.0 : Math::PI
          final = initial
          reverse = wrap_radians(final + Math::PI)
          return GeodesicSolution.new(
            status: :ok,
            converged: true,
            distance: distance,
            initial_forward_azimuth: initial * Ellipsoid::RAD2DEG,
            final_forward_azimuth: final * Ellipsoid::RAD2DEG,
            reverse_azimuth: reverse * Ellipsoid::RAD2DEG,
            iterations: 0,
            distance_error: 0.0,
            longitude_error: 0.0,
            azimuth_error: 0.0,
            message: 'longitude_error is the remaining longitude root error in degrees; azimuth_error is the corresponding angular error in degrees.'
          )
        end

        result = solve_inverse(
          sin_u1,
          cos_u1,
          sin_u2,
          cos_u2,
          lon_diff,
          ellipsoid
        )
        geometry = result[:geometry]
        initial = result[:initial]
        final = result[:final]
        if cos_u1.zero?
          initial = lat1 > 0 ? Math::PI : 0.0
        end
        if cos_u2.zero?
          final = lat2 > 0 ? 0.0 : Math::PI
        end
        reverse = wrap_radians(final + Math::PI)
        azimuth_error = result[:longitude_error]
        distance_error = ellipsoid.semi_minor_axis * result[:longitude_error]
        unless result[:converged]
          return GeodesicSolution.new(
            status: :not_converged,
            converged: false,
            distance: geometry[:distance],
            initial_forward_azimuth: initial * Ellipsoid::RAD2DEG,
            final_forward_azimuth: final * Ellipsoid::RAD2DEG,
            reverse_azimuth: reverse * Ellipsoid::RAD2DEG,
            iterations: result[:iterations],
            distance_error: distance_error,
            longitude_error: result[:longitude_error] * Ellipsoid::RAD2DEG,
            azimuth_error: azimuth_error * Ellipsoid::RAD2DEG,
            message: 'The ellipsoidal iteration did not reach the requested tolerance.'
          )
        end

        GeodesicSolution.new(
          status: :ok,
          converged: true,
          distance: geometry[:distance],
          initial_forward_azimuth: initial * Ellipsoid::RAD2DEG,
          final_forward_azimuth: final * Ellipsoid::RAD2DEG,
          reverse_azimuth: reverse * Ellipsoid::RAD2DEG,
          iterations: result[:iterations],
          distance_error: distance_error,
          longitude_error: result[:longitude_error] * Ellipsoid::RAD2DEG,
          azimuth_error: azimuth_error * Ellipsoid::RAD2DEG,
          message: 'longitude_error is the remaining longitude root error in degrees; azimuth_error is the corresponding angular error in degrees.'
        )
      end

      def direct(point, azimuth, distance, ellipsoid = Ellipsoid::WGS84)
        ellipsoid = validate_ellipsoid(ellipsoid)
        point = validate_point!(point)
        unless finite_number?(azimuth)
          raise ArgumentError, 'azimuth must be a finite number'
        end
        distance = validate_distance!(distance)
        azimuth_radians = wrap_radians(azimuth * Ellipsoid::DEG2RAD)
        lat1 = point.y * Ellipsoid::DEG2RAD
        sin_u1, cos_u1 = reduced_latitude(lat1, ellipsoid.flattening)
        sin_azimuth = Math.sin(azimuth_radians)
        cos_azimuth = Math.cos(azimuth_radians)

        if distance.zero?
          return GeodesicSolution.new(
            status: :same_point,
            converged: true,
            distance: 0.0,
            initial_forward_azimuth: nil,
            final_forward_azimuth: nil,
            reverse_azimuth: nil,
            endpoint: Point.from_x_y(point.x, point.y, point.srid),
            iterations: 0,
            distance_error: 0.0,
            longitude_error: 0.0,
            azimuth_error: 0.0,
            message: 'The distance is zero; azimuth is undefined.'
          )
        end

        sigma1 = Math.atan2(sin_u1, cos_u1 * cos_azimuth)
        sin_alpha_equator = cos_u1 * sin_azimuth
        cos_alpha_equator = Math.sqrt(
          [1.0 - sin_alpha_equator * sin_alpha_equator, 0.0].max
        )
        u_squared = cos_alpha_equator**2 *
          (ellipsoid.semi_major_axis**2 - ellipsoid.semi_minor_axis**2) /
          ellipsoid.semi_minor_axis**2
        a_coefficient = 1.0 + u_squared / 16_384.0 *
          (4096.0 + u_squared * (-768.0 + u_squared * (320.0 - 175.0 * u_squared)))
        b_coefficient = u_squared / 1024.0 *
          (256.0 + u_squared * (-128.0 + u_squared * (74.0 - 47.0 * u_squared)))
        sigma2 = distance / (ellipsoid.semi_minor_axis * a_coefficient) + sigma1
        converged = false
        iterations = 0
        residual = nil
        MAX_ITERATIONS.times do |index|
          iterations = index + 1
          delta = sigma2 - sigma1
          mid = Math.cos(2.0 * sigma1 + delta)
          sin_delta = Math.sin(delta)
          cos_delta = Math.cos(delta)
          correction = b_coefficient * sin_delta *
            (mid + b_coefficient / 4.0 *
              (cos_delta * (-1.0 + 2.0 * mid * mid) -
               b_coefficient / 6.0 * mid * (-3.0 + 4.0 * sin_delta * sin_delta) *
               (-3.0 + 4.0 * mid * mid)))
          next_sigma = distance / (ellipsoid.semi_minor_axis * a_coefficient) +
            sigma1 + correction
          residual = (next_sigma - sigma2).abs
          sigma2 = next_sigma
          if residual <= TOLERANCE
            converged = true
            break
          end
        end

        delta = sigma2 - sigma1
        sin_u2 = cos_alpha_equator * Math.sin(sigma2)
        cos_u2 = Math.sqrt([1.0 - sin_u2 * sin_u2, 0.0].max)
        lat2 = Math.atan2(sin_u2, (1.0 - ellipsoid.flattening) * cos_u2)

        lambda_radians = if cos_u1.zero? || cos_u2.zero?
          0.0
        else
          sin_lambda = sin_azimuth * Math.sin(delta) / cos_u2
          cos_lambda = (Math.cos(delta) - sin_u1 * sin_u2) / (cos_u1 * cos_u2)
          Math.atan2(sin_lambda, cos_lambda)
        end
        root_error = 0.0
        geometry = geometry_at(
          lambda_radians,
          sin_u1,
          cos_u1,
          sin_u2,
          cos_u2,
          ellipsoid
        )
        longitude = point.x + geometry[:longitude] * Ellipsoid::RAD2DEG
        if cos_u1.zero?
          longitude = wrap_degrees(
            point.x + (lat1 > 0 ? azimuth : Math::PI - azimuth) * Ellipsoid::RAD2DEG
          )
        end
        endpoint = Point.from_x_y(
          wrap_degrees(longitude),
          lat2 * Ellipsoid::RAD2DEG,
          point.srid
        )
        final_forward = Math.atan2(
          -cos_u1 * Math.sin(lambda_radians),
          cos_u2 * sin_u1 - sin_u2 * cos_u1 * Math.cos(lambda_radians)
        )
        if cos_u1.zero?
          azimuth_radians = lat1 > 0 ? Math::PI : 0.0
          final_forward = wrap_radians(final_forward)
        end
        final_forward = if lat2.abs == Math::PI / 2.0
          final_forward
        else
          wrap_radians(final_forward)
        end
        reverse = wrap_radians(final_forward + Math::PI)
        distance_error = ellipsoid.semi_minor_axis * residual.to_f
        status = converged && root_error <= TOLERANCE ? :ok : :not_converged
        GeodesicSolution.new(
          status: status,
          converged: status == :ok,
          distance: distance,
          initial_forward_azimuth: wrap_radians(azimuth_radians) * Ellipsoid::RAD2DEG,
          final_forward_azimuth: final_forward * Ellipsoid::RAD2DEG,
          reverse_azimuth: reverse * Ellipsoid::RAD2DEG,
          endpoint: endpoint,
          iterations: iterations,
          distance_error: distance_error,
          longitude_error: root_error * Ellipsoid::RAD2DEG,
          azimuth_error: root_error * Ellipsoid::RAD2DEG,
          message: 'distance_error is the remaining sigma iteration residual in distance units; longitude_error and azimuth_error are zero for the closed longitude expression.'
        )
      end
      class << self
        private

        def validate_ellipsoid(ellipsoid)
          unless ellipsoid.is_a?(Ellipsoid)
            raise ArgumentError, 'ellipsoid must be an Ellipsoid'
          end

          ellipsoid
        end

        def validate_point!(point)
          unless point.is_a?(Point)
            raise ArgumentError, 'point must be a GeoRuby Point'
          end
          unless finite_number?(point.x) && finite_number?(point.y)
            raise ArgumentError, 'point coordinates must be finite numbers'
          end
          unless point.y.between?(-90.0, 90.0)
            raise ArgumentError, 'latitude must be between -90 and 90 degrees'
          end

          point
        end

        def validate_distance!(distance)
          unless distance.is_a?(Numeric) && distance.finite? && distance >= 0
            raise ArgumentError, 'distance must be a finite non-negative number'
          end

          distance.to_f
        end

        def finite_number?(value)
          value.is_a?(Numeric) && value.finite?
        end

        def wrap_degrees(value)
          (value + 180.0) % 360.0 - 180.0
        end

        def wrap_radians(value)
          value % (2.0 * Math::PI)
        end

        def reduced_latitude(latitude, flattening)
          u = Math.atan((1.0 - flattening) * Math.tan(latitude))
          [Math.sin(u), Math.cos(u)]
        end

        def meridian_distance(lat1, lat2, ellipsoid)
          f = ellipsoid.flattening
          n = f / (2.0 - f)
          a_coefficient = 1.0 + n * n / 4.0 + n**4 / 64.0
          b_coefficient = 3.0 / 2.0 * (n - n**3 / 8.0)
          c_coefficient = 15.0 / 16.0 * (n * n - n**4 / 4.0)
          d_coefficient = 35.0 / 48.0 * n**3
          angle = a_coefficient * (lat2 - lat1) -
                  b_coefficient * (Math.sin(2.0 * lat2) - Math.sin(2.0 * lat1)) / 2.0 +
                  c_coefficient * (Math.sin(4.0 * lat2) - Math.sin(4.0 * lat1)) / 4.0 -
                  d_coefficient * (Math.sin(6.0 * lat2) - Math.sin(6.0 * lat1)) / 6.0
          (ellipsoid.semi_major_axis / (1.0 + n)) * angle.abs
        end

        def sphere_point(lambda_radians, sin_u1, cos_u1, sin_u2, cos_u2)
          sin_lambda = Math.sin(lambda_radians)
          cos_lambda = Math.cos(lambda_radians)
          sin_x = cos_u2 * sin_lambda
          cos_x = cos_u1 * sin_u2 - sin_u1 * cos_u2 * cos_lambda
          sin_sigma = Math.hypot(sin_x, cos_x)
          cos_sigma = sin_u1 * sin_u2 + cos_u1 * cos_u2 * cos_lambda
          sigma = Math.atan2(sin_sigma, cos_sigma)
          sin_alpha = sin_sigma.zero? ? 0.0 : cos_u1 * cos_u2 * sin_lambda / sin_sigma
          cos_alpha_squared = 1.0 - sin_alpha * sin_alpha
          cos2_sigma_mid = if cos_alpha_squared < 1e-30
            0.0
          else
            cos_sigma - 2.0 * sin_u1 * sin_u2 / cos_alpha_squared
          end
          cos2_sigma_mid = [[cos2_sigma_mid, -1.0].max, 1.0].min

          {
            sin_x: sin_x,
            cos_x: cos_x,
            sin_sigma: sin_sigma,
            cos_sigma: cos_sigma,
            sigma: sigma,
            sin_alpha: sin_alpha,
            cos_alpha_squared: cos_alpha_squared,
            cos2_sigma_mid: cos2_sigma_mid
          }
        end

        def longitude_at(lambda_radians, geometry, flattening)
          sphere = geometry[:sphere]
          c = flattening / 16.0 * sphere[:cos_alpha_squared] *
              (4.0 + flattening * (4.0 - 3.0 * sphere[:cos_alpha_squared]))
          correction = (1.0 - c) * flattening * sphere[:sin_alpha] *
            (sphere[:sigma] + c * sphere[:sin_sigma] *
              (sphere[:cos2_sigma_mid] + c * sphere[:cos_sigma] *
                (-1.0 + 2.0 * sphere[:cos2_sigma_mid]**2)))
          lambda_radians - correction
        end

        def geometry_at(lambda_radians, sin_u1, cos_u1, sin_u2, cos_u2, ellipsoid)
          sphere = sphere_point(lambda_radians, sin_u1, cos_u1, sin_u2, cos_u2)
          u_squared = sphere[:cos_alpha_squared] *
            (ellipsoid.semi_major_axis**2 - ellipsoid.semi_minor_axis**2) /
            ellipsoid.semi_minor_axis**2
          a_coefficient = 1.0 + u_squared / 16_384.0 *
            (4096.0 + u_squared * (-768.0 + u_squared * (320.0 - 175.0 * u_squared)))
          b_coefficient = u_squared / 1024.0 *
            (256.0 + u_squared * (-128.0 + u_squared * (74.0 - 47.0 * u_squared)))
          mid = sphere[:cos2_sigma_mid]
          delta_sigma = b_coefficient * sphere[:sin_sigma] *
            (mid + b_coefficient / 4.0 *
              (sphere[:cos_sigma] * (-1.0 + 2.0 * mid * mid) -
               b_coefficient / 6.0 * mid * (-3.0 + 4.0 * sphere[:sin_sigma]**2) *
               (-3.0 + 4.0 * mid * mid)))
          distance = ellipsoid.semi_minor_axis * a_coefficient *
            (sphere[:sigma] - delta_sigma)
          longitude = longitude_at(
            lambda_radians,
            { sphere: sphere },
            ellipsoid.flattening
          )
          {
            sphere: sphere,
            distance: distance,
            longitude: longitude,
            a_coefficient: a_coefficient,
            b_coefficient: b_coefficient
          }
        end

        def solve_inverse(sin_u1, cos_u1, sin_u2, cos_u2, longitude_diff, ellipsoid)
          mirror = longitude_diff < 0
          target = mirror ? -longitude_diff : longitude_diff
          lambda_radians = target
          iterations = 0
          converged = false
          previous_residual = nil

          MAX_ITERATIONS.times do |index|
            iterations = index + 1
            sphere = sphere_point(
              lambda_radians,
              sin_u1,
              cos_u1,
              sin_u2,
              cos_u2
            )
            c = ellipsoid.flattening / 16.0 * sphere[:cos_alpha_squared] *
              (4.0 + ellipsoid.flattening *
                (4.0 - 3.0 * sphere[:cos_alpha_squared]))
            next_lambda = target + (1.0 - c) * ellipsoid.flattening *
              sphere[:sin_alpha] *
              (sphere[:sigma] + c * sphere[:sin_sigma] *
                (sphere[:cos2_sigma_mid] + c * sphere[:cos_sigma] *
                  (-1.0 + 2.0 * sphere[:cos2_sigma_mid]**2)))
            residual = (next_lambda - lambda_radians).abs
            lambda_radians = next_lambda
            if residual <= TOLERANCE
              converged = true
              break
            end
            if previous_residual && residual > previous_residual
              lambda_radians = nil
              break
            end
            previous_residual = residual
          end

          if lambda_radians.nil? || !converged
            lower = target
            upper = Math::PI
            MAX_ITERATIONS.times do |index|
              iterations = index + 1
              middle = (lower + upper) / 2.0
              value = geometry_at(
                middle,
                sin_u1,
                cos_u1,
                sin_u2,
                cos_u2,
                ellipsoid
              )[:longitude] - target
              if value < 0
                lower = middle
              else
                upper = middle
              end
            end
            lambda_radians = (lower + upper) / 2.0
            root_error = (upper - lower).abs / 2.0
            converged = root_error <= TOLERANCE
          else
            root_error = (geometry_at(
              lambda_radians,
              sin_u1,
              cos_u1,
              sin_u2,
              cos_u2,
              ellipsoid
            )[:longitude] - target).abs
            converged = root_error <= 1e-10
          end

          geometry = geometry_at(
            lambda_radians,
            sin_u1,
            cos_u1,
            sin_u2,
            cos_u2,
            ellipsoid
          )
          sphere = geometry[:sphere]
          initial = Math.atan2(sphere[:sin_x], sphere[:cos_x])
          final = Math.atan2(
            -cos_u1 * Math.sin(lambda_radians),
            cos_u2 * sin_u1 - sin_u2 * cos_u1 * Math.cos(lambda_radians)
          )
          if mirror
            initial = wrap_radians(2.0 * Math::PI - initial)
            final = wrap_radians(2.0 * Math::PI - final)
          else
            initial = wrap_radians(initial)
            final = wrap_radians(final)
          end

          {
            geometry: geometry,
            initial: initial,
            final: final,
            iterations: iterations,
            converged: converged,
            longitude_error: root_error
          }
        end
      end
    end
  end
end
