# -*- coding: utf-8 -*-
module GeoRuby
  module SimpleFeatures
    # An ellipsoid of revolution used as the reference surface for
    # geodesic computations.  It is defined by its semi-major axis
    # (equatorial radius) +a+ in meters and its flattening +f+.
    class Ellipsoid
      attr_reader :a, :f

      def initialize(a, f)
        unless a.is_a?(Numeric) && a.finite? && a.positive?
          raise ArgumentError,
                "semi-major axis must be a finite positive number, got " \
                  "#{a.inspect}"
        end
        unless f.is_a?(Numeric) && f.finite? && f.abs < 1
          raise ArgumentError,
                "flattening must be a finite number with |f| < 1, got " \
                  "#{f.inspect}"
        end
        @a = a.to_f
        @f = f.to_f
      end

      # Semi-minor (polar) axis, in meters.
      def b
        @a * (1.0 - @f)
      end

      # Square of the second eccentricity: (a^2 - b^2) / b^2.
      def ep2
        f2 = @f * (2.0 - @f)
        f2 / ((1.0 - @f) * (1.0 - @f))
      end
    end

    # The WGS84 ellipsoid, default for the geodesic methods.
    Ellipsoid::WGS84 = Ellipsoid.new(6_378_137.0, 1.0 / 298.257223563)
  end
end

