require 'geo_ruby/simple_features/line_string'

module GeoRuby
  module SimpleFeatures
    # Represents a linear ring, which is a closed line string (see LineString).
    # Currently, no check is performed to verify if the linear ring is really closed.
    class LinearRing < LineString
      SPHERICAL_LATITUDE_EPSILON = 1e-12
      SPHERICAL_DUPLICATE_EPSILON = 1e-15
      SPHERICAL_ANTIPODAL_EPSILON = 1e-12

      def initialize(srid = DEFAULT_SRID, with_z = false, with_m = false)
        super(srid, with_z, with_m)
      end

      # Returns the signed spherical area of this ring.
      #
      # Coordinates are interpreted as longitude/latitude in degrees. Each
      # edge follows the shorter great-circle arc, so rings may cross the
      # antimeridian without wrapping longitude. A counterclockwise ring has
      # positive signed area as viewed from outside the sphere.
      #
      # +radius+ is the sphere radius and determines the area unit. The area
      # is signed: reversing the ring reverses the sign. By default the
      # region not larger than a hemisphere is represented. Pass
      # <tt>larger_region: true</tt> to represent its complementary region.
      # An exact hemisphere boundary necessarily contains antipodal points,
      # so it is rejected by the same uniqueness rule.
      #
      # An open ring is closed cyclically. A repeated final closure point is
      # accepted and ignored, as are other consecutive duplicate points.
      # Rings with fewer than three distinct consecutive vertices and rings
      # contained in one great circle have zero area. Adjacent antipodal
      # points are rejected because infinitely many great circles connect
      # them. Self-intersecting rings are neither detected nor supported;
      # their result is undefined.
      def spherical_area(radius, options = {})
        unless radius.is_a?(Numeric) && radius.real? && radius > 0
          raise ArgumentError, 'radius must be a positive number'
        end

        larger_region = options.fetch(:larger_region, false)
        unless [true, false].include?(larger_region)
          raise ArgumentError, 'larger_region must be true or false'
        end

        vectors = spherical_vectors
        reject_antipodal_edges(vectors)
        return 0.0 if vectors.length < 3

        signed_angle = spherical_signed_angle(vectors)
        if larger_region && !signed_angle.zero?
          signed_angle += signed_angle > 0 ? -4 * Math::PI : 4 * Math::PI
        end

        signed_angle * radius * radius
      end
      alias_method :spherical_signed_area, :spherical_area

      # fix kml export
      alias_method :orig_kml_representation, :kml_representation
      def kml_representation(options = {})
        orig_kml_representation(options).gsub('LineString', 'LinearRing')
      end

      # Does this linear string contain the given point?  We use the
      # algorithm described here:
      #
      # http://www.ecse.rpi.edu/Homepages/wrf/Research/Short_Notes/pnpoly.html
      def contains_point?(point)
        x, y = point.x, point.y
        tuples = @points.zip(@points[1..-1] + [@points[0]])
        crossings =
          tuples.select do |a, b|
            (b.y > y != a.y > y) && (x < (a.x - b.x) * (y - b.y) / (a.y - b.y) + b.x)
          end

        crossings.size % 2 == 1
      end

      private

      def spherical_vectors
        vectors = []

        @points.each do |point|
          unless point.x.is_a?(Numeric) && point.y.is_a?(Numeric) &&
                 point.x.real? && point.y.real?
            raise ArgumentError, 'spherical area requires numeric coordinates'
          end

          latitude = point.y
          if latitude < -90.0 - SPHERICAL_LATITUDE_EPSILON ||
             latitude > 90.0 + SPHERICAL_LATITUDE_EPSILON
            raise ArgumentError, 'latitude must be between -90 and 90 degrees'
          end

          longitude = point.x * Math::PI / 180.0
          latitude = [[-90.0, latitude].max, 90.0].min * Math::PI / 180.0
          cosine_latitude =
            if latitude.abs >= Math::PI / 2.0 - SPHERICAL_DUPLICATE_EPSILON
              0.0
            else
              Math.cos(latitude)
            end
          vector = [
            cosine_latitude * Math.cos(longitude),
            cosine_latitude * Math.sin(longitude),
            Math.sin(latitude)
          ]
          normalize_vector!(vector)

          if vectors.empty? ||
             vector_dot(vectors.last, vector) < 1.0 - SPHERICAL_DUPLICATE_EPSILON
            vectors << vector
          end
        end

        if vectors.length > 1 &&
           vector_dot(vectors.last, vectors.first) >= 1.0 - SPHERICAL_DUPLICATE_EPSILON
          vectors.pop
        end

        vectors
      end

      def reject_antipodal_edges(vectors)
        vectors.each_with_index do |point, index|
          next_point = vectors[(index + 1) % vectors.length]
          next unless next_point

          reject_antipodal_vector_pair(point, next_point)
        end
      end

      def spherical_signed_angle(vectors)
        total_turn = 0.0

        vectors.each_with_index do |current, index|
          previous = vectors[(index - 1) % vectors.length]
          following = vectors[(index + 1) % vectors.length]

          previous_distance = arccosine_clamped(vector_dot(previous, current))
          following_distance = arccosine_clamped(vector_dot(current, following))

          incoming = current.each_with_index.map do |component, coordinate|
            (Math.cos(previous_distance) * component - previous[coordinate]) /
              Math.sin(previous_distance)
          end
          outgoing = following.each_with_index.map do |component, coordinate|
            (component - Math.cos(following_distance) * current[coordinate]) /
              Math.sin(following_distance)
          end
          normalize_vector!(incoming)
          normalize_vector!(outgoing)

          cross_product = vector_cross(incoming, outgoing)
          total_turn += Math.atan2(
            vector_dot(cross_product, current),
            vector_dot(incoming, outgoing)
          )
        end

        signed_angle = 2.0 * Math::PI - total_turn
        four_pi = 4.0 * Math::PI
        two_pi = 2.0 * Math::PI
        while signed_angle > two_pi
          signed_angle -= four_pi
        end
        while signed_angle < -two_pi
          signed_angle += four_pi
        end
        signed_angle
      end

      def arccosine_clamped(value)
        Math.acos([[-1.0, value].max, 1.0].min)
      end

      def reject_antipodal_vector_pair(first, second)
        if vector_dot(first, second) <= -1.0 + SPHERICAL_ANTIPODAL_EPSILON
          raise ArgumentError,
            'adjacent antipodal points do not uniquely define a great-circle arc'
        end
      end

      def vector_dot(first, second)
        first.each_with_index.inject(0.0) do |sum, (component, index)|
          sum + component * second[index]
        end
      end

      def vector_cross(first, second)
        [
          first[1] * second[2] - first[2] * second[1],
          first[2] * second[0] - first[0] * second[2],
          first[0] * second[1] - first[1] * second[0]
        ]
      end

      def normalize_vector!(vector)
        norm = Math.sqrt(
          vector.inject(0.0) { |sum, component| sum + component * component }
        )
        vector.map! { |component| component / norm }
      end
    end
  end
end
