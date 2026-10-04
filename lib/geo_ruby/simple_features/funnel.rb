require 'geo_ruby/simple_features/point'

module GeoRuby
  module SimpleFeatures
    # Computes the shortest polyline ("funnel" path) connecting two
    # points through an ordered corridor of portals.
    #
    # A corridor is described by an ordered list of portals. Each portal
    # is a +[left, right]+ pair of points, oriented as seen when walking
    # from the start of the corridor towards its end. Consecutive
    # portals must share at least one endpoint (they are consecutive
    # edges of a triangle strip) and the left/right orientation must
    # stay consistent along the whole corridor. Zero-width portals
    # (left == right) are allowed and force the path through that point.
    #
    #   path = GeoRuby::SimpleFeatures::Funnel.path(start, goal, portals)
    #
    # +start+, +goal+ and the portal endpoints may be Point instances or
    # +[x, y]+ arrays. The returned path is an array of Points, starting
    # at +start+ and ending at +goal+.
    class Funnel
      attr_reader :portals

      def initialize(portals)
        @portals = portals.map do |portal|
          left, right = portal
          [Funnel.to_point(left), Funnel.to_point(right)]
        end
        validate!
      end

      # Builds a funnel for +portals+ and returns the shortest path
      # from +start+ to +goal+ inside the corridor.
      def self.path(start, goal, portals)
        new(portals).path(start, goal)
      end

      # Shortest polyline from +start+ to +goal+ staying inside the
      # corridor. Returns an array of Points.
      def path(start, goal)
        start = Funnel.to_point(start)
        goal = Funnel.to_point(goal)
        return [start, goal] if portals.empty?

        lefts = portals.map(&:first)
        rights = portals.map(&:last)
        # The goal acts as a final zero-width portal.
        lefts << goal
        rights << goal

        output = [start]
        apex = start
        portal_left = lefts[0]
        portal_right = rights[0]
        left_index = 0
        right_index = 0

        i = 1
        while i < lefts.length
          # Tighten the right side of the funnel.
          point = rights[i]
          if Funnel.triarea2(apex, portal_right, point) >= 0.0
            if Funnel.same_point?(apex, portal_right) ||
               Funnel.triarea2(apex, portal_left, point) < 0.0
              portal_right = point
              right_index = i
            else
              # The right side crossed over the left side: the path is
              # taut at the left point, which becomes the new apex.
              apex = portal_left
              output << apex unless Funnel.same_point?(output.last, apex)
              portal_left = apex
              portal_right = apex
              i = left_index + 1
              left_index = i - 1
              right_index = i - 1
              next
            end
          end

          # Tighten the left side of the funnel.
          point = lefts[i]
          if Funnel.triarea2(apex, portal_left, point) <= 0.0
            if Funnel.same_point?(apex, portal_left) ||
               Funnel.triarea2(apex, portal_right, point) > 0.0
              portal_left = point
              left_index = i
            else
              # The left side crossed over the right side.
              apex = portal_right
              output << apex unless Funnel.same_point?(output.last, apex)
              portal_left = apex
              portal_right = apex
              i = right_index + 1
              left_index = i - 1
              right_index = i - 1
              next
            end
          end

          i += 1
        end

        output << goal unless Funnel.same_point?(output.last, goal)
        output
      end

      # Twice the signed area of triangle (a, b, c); positive when +c+
      # lies to the left of the directed line a -> b.
      def self.triarea2(a, b, c)
        (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x)
      end

      def self.same_point?(a, b)
        a.x == b.x && a.y == b.y
      end

      def self.to_point(point)
        case point
        when Point
          point
        when Array
          Point.from_x_y(point[0], point[1])
        else
          if point.respond_to?(:x) && point.respond_to?(:y)
            Point.from_x_y(point.x, point.y)
          else
            raise ArgumentError, 'cannot use given object as a 2D point'
          end
        end
      end

      private

      # Checks the continuity of the corridor: consecutive portals must
      # share at least one endpoint, and the shared endpoint must stay
      # on the same side of the corridor (the left/right orientation may
      # not flip from one portal to the next). Zero-width portals make
      # the side ambiguous and are exempt from the orientation check.
      def validate!
        portals.each_cons(2) do |(left_a, right_a), (left_b, right_b)|
          same_side = Funnel.same_point?(left_a, left_b) ||
                      Funnel.same_point?(right_a, right_b)
          flipped = Funnel.same_point?(left_a, right_b) ||
                    Funnel.same_point?(right_a, left_b)
          unless same_side || flipped
            raise ArgumentError,
                  'portals must be continuous: ' \
                  'consecutive portals must share an endpoint'
          end
          if flipped && !same_side &&
             !Funnel.same_point?(left_a, right_a) &&
             !Funnel.same_point?(left_b, right_b)
            raise ArgumentError,
                  'inconsistent portal orientation: ' \
                  'shared endpoint switches side'
          end
        end
      end
    end
  end
end
