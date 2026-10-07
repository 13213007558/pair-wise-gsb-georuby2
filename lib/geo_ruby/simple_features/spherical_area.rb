require 'geo_ruby/simple_features/linear_ring'
require 'geo_ruby/simple_features/polygon'

module GeoRuby
  module SimpleFeatures
    # Spherical (great-circle) area for geographic rings and polygons,
    # computed in pure Ruby -- no GIS library is required.
    #
    # Coordinates are interpreted as longitude (x) / latitude (y) in
    # degrees. Every edge is interpreted as the *shorter* great-circle
    # arc between its endpoints, and the sphere radius is supplied by
    # the caller (the result is expressed in radius^2 units, e.g.
    # square meters when the radius is given in meters).
    #
    # Semantics and supported input scope:
    #
    # * Direction matters: the signed area is positive when the ring is
    #   traversed counter-clockwise as seen from outside the sphere
    #   (the enclosed region lies to the left of the boundary) and
    #   negative otherwise. Reversing a ring negates the result.
    # * Antimeridian: computation is done on 3D unit vectors, never on
    #   longitude differences, so rings crossing the 180th meridian
    #   follow the short arc instead of jumping the wrong way around
    #   the globe.
    # * Antipodal vertices: two consecutive vertices that are (near-)
    #   antipodal do not define a unique shorter arc, so an
    #   ArgumentError is raised instead of silently picking an
    #   arbitrary great-circle plane. The same applies to a vertex
    #   antipodal to the fan origin used internally.
    # * Regions larger than a hemisphere: a ring bounds two
    #   complementary regions, so when the magnitude of the result
    #   exceeds a hemisphere the caller must pass an explicit +rule+
    #   (+:signed+, +:left+ or +:smaller+); without one an
    #   ArgumentError is raised.
    # * A repeated closing point (first == last) is tolerated and
    #   ignored; it does not change the result.
    # * Rings with fewer than 3 distinct vertices, or collapsed onto a
    #   single great circle, are zero-area rings and return 0.0.
    # * Self-intersecting rings are not validated or rejected; the
    #   result is the algebraic (signed) sum of the lobes, so lobes
    #   wound in opposite directions cancel out.
    # * Polygon holes: the polygon total is the plain sum of the
    #   signed ring areas. Holes are expected to be wound opposite to
    #   the outer ring so they subtract; the signed values are never
    #   absolutized, so a hole wound like the outer ring adds its
    #   area instead of subtracting it.
    module SphericalArea
      # Dot product below which two unit vectors are considered
      # (near-)antipodal and the arc between them is rejected.
      ANTIPODAL_DOT_THRESHOLD = -1.0 + 1e-10

      # Rules accepted by #spherical_area for resolving how the signed
      # result is reported, mandatory once a region exceeds a hemisphere.
      RULES = %i[signed left smaller].freeze

      # Area of the geometry on a sphere of the given +radius+.
      #
      # Without a +rule+, returns the signed area and raises
      # ArgumentError when its magnitude exceeds a hemisphere. Rules:
      #
      # * +:signed+   -- signed area in (-4*pi*R^2, 4*pi*R^2); the sign
      #   follows the traversal direction.
      # * +:left+     -- area of the region to the left of the
      #   traversal direction, in [0, 4*pi*R^2).
      # * +:smaller+  -- area of the smaller of the two complementary
      #   regions bounded by the ring, in [0, 2*pi*R^2].
      def spherical_area(radius, rule: nil)
        SphericalArea.apply_rule(spherical_signed_area(radius), radius, rule)
      end

      # Applies a reporting rule to a signed spherical area. See
      # #spherical_area for the meaning of the rules.
      def self.apply_rule(signed_area, radius, rule)
        full_sphere = 4.0 * Math::PI * radius * radius
        if rule.nil?
          hemisphere = full_sphere / 2.0
          if signed_area.abs > hemisphere
            raise ArgumentError,
                  'area exceeds a hemisphere; pass an explicit rule (:signed, :left or :smaller)'
          end
          signed_area
        else
          case rule
          when :signed then signed_area
          when :left then signed_area.negative? ? signed_area + full_sphere : signed_area
          when :smaller then [signed_area.abs, full_sphere - signed_area.abs].min
          else
            raise ArgumentError,
                  "unknown spherical area rule #{rule.inspect}; expected one of #{RULES.inspect}"
          end
        end
      end

      private

      # Unit vector on the sphere for a lon/lat point (degrees).
      def spherical_unit_vector(point)
        lon = point.x * Math::PI / 180.0
        lat = point.y * Math::PI / 180.0
        cos_lat = Math.cos(lat)
        [cos_lat * Math.cos(lon), cos_lat * Math.sin(lon), Math.sin(lat)]
      end

      def spherical_dot(a, b)
        a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
      end

      def spherical_check_not_antipodal!(a, b)
        return unless spherical_dot(a, b) < ANTIPODAL_DOT_THRESHOLD

        raise ArgumentError,
              '(near-)antipodal vertices: the shorter great-circle arc is not unique'
      end

      # Signed area (steradians) of the spherical triangle a-b-c given
      # as unit vectors, using the Van Oosterom-Strackee formula.
      def spherical_triangle_signed_steradians(a, b, c)
        triple = a[0] * (b[1] * c[2] - b[2] * c[1]) -
                 a[1] * (b[0] * c[2] - b[2] * c[0]) +
                 a[2] * (b[0] * c[1] - b[1] * c[0])
        denominator = 1.0 + spherical_dot(a, b) + spherical_dot(b, c) +
                      spherical_dot(a, c)
        2.0 * Math.atan2(triple, denominator)
      end

      # Signed steradians of a ring given as an array of unit vectors
      # (unclosed), via fan triangulation from the first vertex.
      def spherical_ring_signed_steradians(vectors)
        return 0.0 if vectors.length < 3

        origin = vectors.first
        vectors.each_cons(2) { |a, b| spherical_check_not_antipodal!(a, b) }
        spherical_check_not_antipodal!(vectors.last, origin)
        sum = 0.0
        (1...(vectors.length - 1)).each do |i|
          spherical_check_not_antipodal!(origin, vectors[i])
          sum += spherical_triangle_signed_steradians(origin, vectors[i], vectors[i + 1])
        end
        sum
      end
    end

    class LinearRing
      include SphericalArea

      # Signed spherical area of the ring on a sphere of the given
      # +radius+, in (-4*pi*R^2, 4*pi*R^2). Positive when the enclosed
      # region lies to the left of the traversal direction. See
      # SphericalArea for the supported input scope.
      def spherical_signed_area(radius)
        pts = points
        pts = pts[0...-1] if pts.length > 1 && pts.first == pts.last
        vectors = pts.map { |point| spherical_unit_vector(point) }
        spherical_ring_signed_steradians(vectors) * radius * radius
      end
    end

    class Polygon
      include SphericalArea

      # Signed spherical area of the polygon on a sphere of the given
      # +radius+: the plain sum of the signed areas of the outer ring
      # and the holes. Holes wound opposite to the outer ring subtract
      # from the total; no absolute values are taken. See SphericalArea
      # for the supported input scope.
      def spherical_signed_area(radius)
        rings.sum { |ring| ring.spherical_signed_area(radius) }
      end
    end
  end
end
