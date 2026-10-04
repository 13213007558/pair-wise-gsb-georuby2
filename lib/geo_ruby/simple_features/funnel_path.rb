require 'geo_ruby/simple_features/line_string'
require 'geo_ruby/simple_features/point'

module GeoRuby
  module SimpleFeatures
    # Shortest polyline through an ordered sequence of portal segments.
    class FunnelPath
      class Error < ArgumentError; end

      EPSILON = 1e-12

      attr_reader :start_point, :end_point, :portals

      def self.shortest_path(start_point, end_point, portals)
        new(start_point, end_point, portals).path
      end
      class << self
        alias_method :generate, :shortest_path
        alias_method :call, :shortest_path
      end

      def initialize(start_point, end_point, portals)
        @start_point = coerce_point(start_point, DEFAULT_SRID)
        @end_point = coerce_point(end_point, @start_point.srid)
        @portals = normalize_portals(portals)
      end

      def path
        LineString.from_points(shortest_points, start_point.srid)
      end
      alias_method :to_line_string, :path

      private

      def shortest_points
        return [start_point] if same_point?(start_point, end_point)

        nodes = [{ point: start_point, index: -1 }]
        portals.each_with_index do |(left, right), index|
          nodes << { point: left, index: index }
          nodes << { point: right, index: index } unless same_point?(left, right)
        end
        nodes << { point: end_point, index: portals.length }

        distances = Array.new(nodes.length) { Float::INFINITY }
        previous = Array.new(nodes.length)
        distances[0] = 0.0

        nodes.each_with_index do |source, source_index|
          next if distances[source_index].infinite?

          nodes.each_with_index do |target, target_index|
            next unless target_index > source_index
            next unless target[:index] > source[:index]
            next unless visible?(source, target)

            distance = distances[source_index] +
                       source[:point].euclidian_distance(target[:point])
            next unless distance < distances[target_index]

            distances[target_index] = distance
            previous[target_index] = source_index
          end
        end

        if distances.last.infinite?
          raise Error, 'no path through the ordered portals exists'
        end

        reconstruct_path(nodes, previous)
      end

      def visible?(source, target)
        (source[:index] + 1...target[:index]).all? do |portal_index|
          left, right = portals[portal_index]
          segments_intersect?(source[:point], target[:point], left, right)
        end
      end

      def reconstruct_path(nodes, previous)
        path = []
        node_index = nodes.length - 1

        while node_index
          point = nodes[node_index][:point]
          path.unshift(point) if path.empty? || !same_point?(path.first, point)
          node_index = previous[node_index]
        end

        path
      end

      def normalize_portals(portal_list)
        unless portal_list.respond_to?(:each)
          raise Error, 'portals must respond to #each'
        end

        portals = portal_list.map.with_index { |portal, index| coerce_portal(portal, index) }
        validate_continuity(portals)
        portals
      end

      def coerce_portal(portal, index)
        endpoints =
          if portal.is_a?(LineString)
            portal.points
          elsif portal.respond_to?(:to_a)
            portal.to_a
          else
            raise Error, "portal ##{index} must contain two endpoints"
          end

        unless endpoints.length == 2
          raise Error, "portal ##{index} must contain two endpoints"
        end

        endpoints.map { |point| coerce_point(point, start_point.srid) }
      end

      def validate_continuity(portals)
        portals.each_cons(2).with_index do |(previous, current), index|
          unless portals_share_endpoint?(previous, current)
            raise Error, "portals ##{index} and ##{index + 1} do not share an endpoint"
          end
        end
      end

      def coerce_point(point, srid)
        case point
        when Point
          Point.from_x_y(Float(point.x), Float(point.y), point.srid)
        when Array
          Point.from_x_y(Float(point[0]), Float(point[1]), srid)
        else
          raise Error, 'expected a Point or [x, y] coordinates'
        end
      rescue ArgumentError, TypeError, NoMethodError
        raise Error, 'endpoints must contain finite x/y coordinates'
      end

      def portals_share_endpoint?(first, second)
        first.any? { |point| second.any? { |other| same_point?(point, other) } }
      end

      def same_point?(first, second)
        (first.x - second.x).abs <= epsilon(first, second) &&
          (first.y - second.y).abs <= epsilon(first, second)
      end

      def segments_intersect?(first, second, third, fourth)
        orientations = [
          cross(third, fourth, first),
          cross(third, fourth, second),
          cross(first, second, third),
          cross(first, second, fourth)
        ]

        return true if point_on_segment?(third, fourth, first)
        return true if point_on_segment?(third, fourth, second)
        return true if point_on_segment?(first, second, third)
        return true if point_on_segment?(first, second, fourth)

        opposite_sides?(orientations[0], orientations[1]) &&
          opposite_sides?(orientations[2], orientations[3])
      end

      def opposite_sides?(first, second)
        (first.negative? && second.positive?) ||
          (first.positive? && second.negative?)
      end

      def point_on_segment?(first, second, point)
        return false unless near_zero?(cross(first, second, point),
                                    epsilon(first, second, point))

        point.x >= [first.x, second.x].min - epsilon(first, second) &&
          point.x <= [first.x, second.x].max + epsilon(first, second) &&
          point.y >= [first.y, second.y].min - epsilon(first, second) &&
          point.y <= [first.y, second.y].max + epsilon(first, second)
      end

      def cross(first, second, third)
        (second.x - first.x) * (third.y - first.y) -
          (second.y - first.y) * (third.x - first.x)
      end

      def epsilon(*points)
        scale = points.map { |point| point.x.abs + point.y.abs }.max
        EPSILON * [1.0, scale].max
      end

      def near_zero?(value, tolerance)
        value.abs <= tolerance
      end
    end

    class LineString < Geometry
      def self.funnel_path(start_point, end_point, portals)
        FunnelPath.shortest_path(start_point, end_point, portals)
      end
    end
  end
end
