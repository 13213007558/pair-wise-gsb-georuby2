require 'set'

module GeoRuby
  module SimpleFeatures
    class PolygonSetSimplifier
      Result = Struct.new(:multi_polygon, :warnings, keyword_init: true)

      Ring = Struct.new(:polygon_index, :ring_index, :points, :closed,
                        :srid, :with_z, :with_m, keyword_init: true)
      Edge = Struct.new(:id, :polygon_index, :ring_index, :start_key, :end_key,
                        :start_point, :end_point, :shared, :candidate,
                        keyword_init: true)
      SharedGroup = Struct.new(:id, :edges, :start_key, :end_key, keyword_init: true)
      SharedArc = Struct.new(:nodes, :occurrences, keyword_init: true)

      COLLINEAR_EPSILON = 1e-12

      def initialize(polygons, epsilon = 1, anchors: [])
        @polygons = Array(polygons)
        @epsilon = epsilon.to_f
        @anchors = Array(anchors).map { |anchor| normalize_anchor(anchor) }
        @warnings = []
      end

      def simplify
        build_rings
        build_ring_point_maps
        validate_anchors
        build_edges
        identify_shared_boundaries
        build_ring_representatives
        simplify_shared_arcs
        simplify_nonshared_rings
        build_result
      end

      def self.simplify(polygons, epsilon = 1, anchors: [])
        new(polygons, epsilon, anchors: anchors).simplify
      end

      private

      def normalize_anchor(anchor)
        case anchor
        when Point then [anchor.x, anchor.y]
        when Array then [anchor[0].to_f, anchor[1].to_f]
        else raise ArgumentError, 'anchors must be Point objects or [x, y] arrays'
        end
      end

      def build_rings
        @rings = []
        @polygons.each_with_index do |polygon, polygon_index|
          polygon.rings.each_with_index do |ring, ring_index|
            points = ring.points.dup
            closed = ring.is_a?(LinearRing) || points.first == points.last
            points.pop if closed && points.length > 1 && points.first == points.last
            @rings << Ring.new(
              polygon_index: polygon_index,
              ring_index: ring_index,
              points: points,
              closed: closed,
              srid: polygon.srid,
              with_z: polygon.with_z?,
              with_m: polygon.with_m?
            )
          end
        end
      end

      def build_ring_point_maps
        @ring_point_maps = Hash.new { |hash, key| hash[key] = {} }
        @rings.each do |ring|
          ring.points.each do |point|
            @ring_point_maps[ring_key(ring)][point_key(point)] = point
          end
        end
      end

      def validate_anchors
        @anchors.each do |anchor|
          on_boundary = @ring_point_maps.values.any? { |points| points.key?(anchor) } ||
                        @rings.any? do |ring|
                          ring.points.each_cons(2).any? { |a, b| point_on_segment?(anchor, a, b) }
                        end
          raise ArgumentError, "anchor #{anchor.inspect} is not on a polygon boundary" unless on_boundary
        end
      end

      def build_edges
        @edges = []
        @ring_edges = Hash.new { |hash, key| hash[key] = [] }
        @node_graph_incidences = Hash.new(0)
        @node_polygons = Hash.new { |hash, key| hash[key] = Set.new }
        @all_points = {}
        @next_edge_id = -1

        @rings.each { |ring| ring.points.each { |point| @all_points[point_key(point)] = canonical_point(point) } }
        @anchors.each { |key| @all_points[key] = Point.from_x_y(key[0], key[1]) }

        edge_id = 0
        @rings.each do |ring|
          boundary_edges(ring).each do |start_point, end_point|
            split_segment_keys(start_point, end_point).each_cons(2) do |start_key, end_key|
              store_ring_point(ring, start_key, start_point, end_point)
              store_ring_point(ring, end_key, start_point, end_point)
              edge = Edge.new(
                id: edge_id,
                polygon_index: ring.polygon_index,
                ring_index: ring.ring_index,
                start_key: start_key,
                end_key: end_key,
                start_point: @all_points[start_key],
                end_point: @all_points[end_key],
                shared: false,
                candidate: false
              )
              @edges << edge
              @ring_edges[ring_key(ring)] << edge
              [start_key, end_key].each do |key|
                @node_polygons[key] << ring.polygon_index
              end
              edge_id += 1
            end
          end
        end
        @edges.group_by { |edge| geometry_key(edge.start_key, edge.end_key) }.each_value do |edges|
          edge = edges.first
          @node_graph_incidences[edge.start_key] += 1
          @node_graph_incidences[edge.end_key] += 1
        end
        @next_edge_id = edge_id - 1
        @retired_edge_ids = Set.new
      end

      def boundary_edges(ring)
        edges = ring.points.each_cons(2).to_a
        edges << [ring.points.last, ring.points.first] if ring.closed && ring.points.length > 2
        edges
      end

      def split_segment_keys(start_point, end_point)
        start = point_key(start_point)
        finish = point_key(end_point)
        candidates = @all_points.keys + @anchors
        keys = candidates.uniq.select do |key|
          key != start && key != finish && exactly_on_segment?(key, start_point, end_point)
        end
        [start] + keys.sort_by { |key| segment_ratio(start_point, end_point, key) } + [finish]
      end

      def exactly_on_segment?(key, start_point, end_point)
        x, y = key
        dx = end_point.x - start_point.x
        dy = end_point.y - start_point.y
        return false if dx.zero? && dy.zero?

        cross = (x - start_point.x) * dy - (y - start_point.y) * dx
        cross.zero? && raw_segment_ratio(start_point, end_point, key).between?(0.0, 1.0)
      end

      def collinear_with_segment?(key, start_point, end_point)
        x, y = key
        dx = end_point.x - start_point.x
        dy = end_point.y - start_point.y
        length_squared = dx * dx + dy * dy
        return false if length_squared.zero?

        cross = ((x - start_point.x) * dy - (y - start_point.y) * dx).abs
        return false unless cross <= COLLINEAR_EPSILON * [1.0, length_squared].max

        raw_segment_ratio(start_point, end_point, key).between?(-COLLINEAR_EPSILON, 1 + COLLINEAR_EPSILON)
      end

      def segment_ratio(start_point, end_point, key)
        (raw_segment_ratio(start_point, end_point, key)).clamp(0.0, 1.0)
      end

      def raw_segment_ratio(start_point, end_point, key)
        dx = end_point.x - start_point.x
        dy = end_point.y - start_point.y
        length_squared = dx * dx + dy * dy
        return 0.0 if length_squared.zero?

        dot = (key[0] - start_point.x) * dx + (key[1] - start_point.y) * dy
        dot / length_squared
      end

      def point_on_segment?(key, start_point, end_point)
        return false unless collinear_with_segment?(key, start_point, end_point)

        ratio = segment_ratio(start_point, end_point, key)
        ratio > COLLINEAR_EPSILON && ratio < 1 - COLLINEAR_EPSILON
      end

      def store_ring_point(ring, key, start_point, end_point)
        map = @ring_point_maps[ring_key(ring)]
        return if map.key?(key)

        ratio = segment_ratio(start_point, end_point, key)
        point =
          if ring.with_z && ring.with_m
            Point.from_x_y_z_m(key[0], key[1],
                               interpolate(start_point.z, end_point.z, ratio),
                               interpolate(start_point.m, end_point.m, ratio), ring.srid)
          elsif ring.with_z
            Point.from_x_y_z(key[0], key[1], interpolate(start_point.z, end_point.z, ratio), ring.srid)
          elsif ring.with_m
            Point.from_x_y_m(key[0], key[1], interpolate(start_point.m, end_point.m, ratio), ring.srid)
          else
            Point.from_x_y(key[0], key[1], ring.srid)
          end
        map[key] = point
      end

      def interpolate(start_value, end_value, ratio)
        start_value + (end_value - start_value) * ratio
      end

      def identify_shared_boundaries
        @shared_groups = []
        grouped = @edges.group_by { |edge| geometry_key(edge.start_key, edge.end_key) }
        grouped.each do |_key, edges|
          next unless edges.length > 1

          start_key, end_key = geometry_key(edges.first.start_key, edges.first.end_key)
          group = SharedGroup.new(id: @shared_groups.length, edges: edges,
                                  start_key: start_key, end_key: end_key)
          edges.each { |edge| edge.shared = true }
          @shared_groups << group
        end
        build_shared_arcs
      end

      def build_shared_arcs
        adjacency = Hash.new { |hash, key| hash[key] = [] }
        @shared_groups.each do |group|
          adjacency[group.start_key] << group
          adjacency[group.end_key] << group
        end

        unvisited = @shared_groups.to_set
        @shared_arcs = []
        while unvisited.any?
          start_group = unvisited.first
          endpoint = [start_group.start_key, start_group.end_key].find do |key|
            degree = adjacency[key].length
            degree == 1 || degree != 2
          end
          start_key = endpoint || start_group.start_key
          group = start_group
          current_key = start_key
          nodes = [start_key]
          occurrences = Hash.new { |hash, key| hash[key] = { edges: [] } }

          loop do
            unvisited.delete(group)
            next_key = group.start_key == current_key ? group.end_key : group.start_key
            group.edges.each do |edge|
              occurrences[ring_key_for_edge(edge)][:edges] << edge
            end
            nodes << next_key
            connected = adjacency[next_key].reject { |item| item == group || !unvisited.include?(item) }
            break if connected.length != 1

            group = connected.first
            current_key = next_key
          end

          occurrences.each_value do |occurrence|
            occurrence[:edges].sort_by!(&:id)
            occurrence[:nodes] = edge_chain_nodes(occurrence[:edges])
          end
          @shared_arcs << SharedArc.new(nodes: nodes, occurrences: occurrences)
        end
      end

      def edge_chain_nodes(edges)
        return [] if edges.empty?

        remaining = edges.dup
        current = remaining.shift
        nodes = [current.start_key, current.end_key]
        until remaining.empty?
          edge = remaining.find { |item| item.start_key == nodes.last } ||
                 remaining.find { |item| item.end_key == nodes.last }
          break unless edge

          remaining.delete(edge)
          nodes << (edge.start_key == nodes.last ? edge.end_key : edge.start_key)
        end
        nodes
      end

      def point_key(point)
        [point.x, point.y]
      end

      def canonical_point(point)
        Point.from_x_y(point.x, point.y, point.srid)
      end

      def ring_key(ring)
        [ring.polygon_index, ring.ring_index]
      end

      def ring_key_for_edge(edge)
        [edge.polygon_index, edge.ring_index]
      end

      def geometry_key(start_key, end_key)
        [[start_key, end_key], [end_key, start_key]].min
      end

      def build_ring_representatives
        original_cycles = {}
        @ring_edges.each do |key, edges|
          original_cycles[key] = cycle_keys(edges)
        end

        @ring_samples = {}
        @original_contains = {}
        @ring_edges.each do |ring_key, edges|
          @ring_samples[ring_key] = interior_edge_samples(cycle_keys(edges))
        end
        @ring_edges.each do |ring_key, edges|
          @original_contains[ring_key] = {}
          cycle = cycle_keys(edges)
          @ring_samples.each_value do |samples|
            samples.each do |sample|
              @original_contains[ring_key][sample] = point_in_ring?(sample, cycle)
            end
          end
        end
      end

      def interior_edge_samples(cycle)
        coordinates = cycle[0...-1]
        area = ring_signed_area(cycle)
        samples = coordinates.each_with_index.map do |point, index|
          previous = coordinates[(index - 1) % coordinates.length]
          next_point = coordinates[(index + 1) % coordinates.length]
          local_area = signed_area(previous, point, next_point)
          next unless area.positive? ? local_area.negative? : local_area.positive?

          centroid = [(previous[0] + point[0] + next_point[0]) / 3.0,
                      (previous[1] + point[1] + next_point[1]) / 3.0]
          inner_midpoint = [(previous[0] + next_point[0]) / 2.0,
                            (previous[1] + next_point[1]) / 2.0]
          candidate_samples = [centroid]
          4.times do |index|
            ratio = (index + 1).to_f / 5.0
            candidate_samples << [
              point[0] + (inner_midpoint[0] - point[0]) * ratio,
              point[1] + (inner_midpoint[1] - point[1]) * ratio
            ]
          end
          candidate_samples.find { |sample| point_in_ring?(sample, cycle) }
        end
          .compact
        samples << interior_seed(cycle)
        samples.compact.uniq
      end

      def ring_signed_area(cycle)
        cycle.each_cons(2).sum do |start_key, end_key|
          start_key[0] * end_key[1] - end_key[0] * start_key[1]
        end / 2.0
      end

      def signed_area(first, second, third)
        (second[0] - first[0]) * (third[1] - first[1]) -
          (second[1] - first[1]) * (third[0] - first[0])
      end

      def interior_seed(cycle)
        coordinates = cycle[0...-1]
        min_x, max_x = coordinates.map(&:first).minmax
        min_y, max_y = coordinates.map(&:last).minmax
        (0..20).each do |attempt|
          x = min_x + (max_x - min_x) * ((attempt * 37) % 101) / 100.0
          y = min_y + (max_y - min_y) * ((attempt * 53) % 101) / 100.0
          return [x, y] if point_in_ring?([x, y], cycle)
        end
        nil
      end

      def cycle_keys(edges)
        return [] if edges.empty?

        remaining = edges.dup
        first = remaining.shift
        keys = [first.start_key, first.end_key]
        until remaining.empty?
          edge = remaining.find { |item| item.start_key == keys.last } ||
                 remaining.find { |item| item.end_key == keys.last }
          break unless edge

          remaining.delete(edge)
          keys << (edge.start_key == keys.last ? edge.end_key : edge.start_key)
        end
        keys
      end

      def point_in_ring?(key, cycle)
        inside = false
        x, y = key
        cycle.each_cons(2) do |start_key, end_key|
          sx, sy = start_key
          ex, ey = end_key
          next if point_on_segment?(key, Point.from_x_y(sx, sy), Point.from_x_y(ex, ey))

          if ((sy <= y && ey > y) || (ey <= y && sy > y)) &&
             x < (sx + (ex - sx) * (y - sy) / (ey - sy))
            inside = !inside
          end
        end
        inside
      end

      def distance_to_segment(key, start_key, end_key)
        point = Point.from_x_y(key[0], key[1])
        point.orthogonal_distance(Point.from_x_y(start_key[0], start_key[1]),
                                  Point.from_x_y(end_key[0], end_key[1]))
      end

      def simplify_shared_arcs
        @shared_arcs.each do |arc|
          mandatory = (0...arc.nodes.length).select do |index|
            key = arc.nodes[index]
            index.zero? || index == arc.nodes.length - 1 ||
              @anchors.include?(key) || junction_node?(key)
          end
          mandatory.each_cons(2) do |start_index, end_index|
            simplify_shared_interval(arc, start_index, end_index)
          end
        end
      end

      def junction_node?(key)
        @anchors.include?(key) || @node_graph_incidences[key] != 2 || @node_polygons[key].length > 2
      end

      def simplify_shared_interval(arc, start_index, end_index)
        return if end_index - start_index <= 1

        farthest_index = (start_index + 1...end_index).max_by do |index|
          distance_to_segment(arc.nodes[index], arc.nodes[start_index], arc.nodes[end_index])
        end
        farthest_distance = distance_to_segment(arc.nodes[farthest_index],
                                                arc.nodes[start_index], arc.nodes[end_index])

        if farthest_distance <= @epsilon
          proposals = shared_interval_proposals(arc, start_index, end_index)
          if proposals.any? && valid_proposals?(proposals)
            apply_proposals(proposals)
          else
            retain_interval(arc.nodes[start_index], arc.nodes[end_index], farthest_distance)
            simplify_shared_interval(arc, start_index, farthest_index)
            simplify_shared_interval(arc, farthest_index, end_index)
          end
        else
          simplify_shared_interval(arc, start_index, farthest_index)
          simplify_shared_interval(arc, farthest_index, end_index)
        end
      end

      def shared_interval_proposals(arc, start_index, end_index)
        proposals = Hash.new { |hash, key| hash[key] = { retired: Set.new, added: [] } }

        arc.occurrences.each do |ring_key, occurrence|
          global_start = arc.nodes[start_index]
          global_end = arc.nodes[end_index]
          occurrence_start = occurrence[:nodes].index(global_start)
          occurrence_end = occurrence[:nodes].index(global_end)
          next unless occurrence_start && occurrence_end

          start_position, end_position = [occurrence_start, occurrence_end].sort
          (start_position...end_position).each do |edge_index|
            edge = occurrence[:edges][edge_index]
            proposals[ring_key][:retired] << edge.id if edge
          end
          ordered_nodes = occurrence_start <= occurrence_end ?
                             [global_start, global_end] : [global_end, global_start]
          edge = Edge.new(
            id: nil,
            polygon_index: ring_key[0],
            ring_index: ring_key[1],
            start_key: ordered_nodes.first,
            end_key: ordered_nodes.last,
            start_point: @all_points[ordered_nodes.first],
            end_point: @all_points[ordered_nodes.last],
            shared: true,
            candidate: true
          )
          proposals[ring_key][:added] << edge
        end
        proposals
      end

      def retain_interval(start_key, end_key, distance)
        @warnings << {
          reason: 'topology_constraint',
          coordinates: [start_key, end_key],
          tolerance: @epsilon,
          deviation: distance,
          message: 'kept interval because the simplified chord would cross a boundary, touch a separate face, or change hole ownership'
        }
      end

      def valid_proposals?(proposals)
        temporary_edges = {}
        affected_rings = proposals.keys

        @ring_edges.each do |ring_key, edges|
          proposal = proposals[ring_key]
          temporary_edges[ring_key] = if proposal
                                       edges.reject { |edge| @retired_edge_ids.include?(edge.id) || proposal[:retired].include?(edge.id) } + proposal[:added]
                                     else
                                       edges.reject { |edge| @retired_edge_ids.include?(edge.id) }
                                     end
        end

        return false unless containment_preserved?(affected_rings, temporary_edges)
        return false unless intersections_preserved?(proposals, temporary_edges)

        true
      end

      def containment_preserved?(affected_rings, temporary_edges)
        affected_rings.all? do |ring_key|
          cycle = cycle_keys(temporary_edges[ring_key])
          next false if cycle.length < 4 || cycle.first != cycle.last

          @ring_samples.all? do |_sample_ring_key, samples|
            samples.all? do |sample|
              point_in_ring?(sample, cycle) == @original_contains[ring_key][sample]
            end
          end
        end
      end

      def intersections_preserved?(proposals, temporary_edges)
        added = proposals.values.flat_map { |proposal| proposal[:added] }
        added.each do |edge|
          temporary_edges.each do |ring_key, edges|
            edges.each do |other|
              next if edge.equal?(other)
              relation = segment_relation(edge, other)
              next if relation == :none
              if relation == :endpoint_touch && edge.polygon_index != other.polygon_index
                shared_endpoint = ([edge.start_key, edge.end_key] &
                                  [other.start_key, other.end_key]).any? do |key|
                  @node_polygons[key].include?(edge.polygon_index) &&
                    @node_polygons[key].include?(other.polygon_index)
                end
                next if shared_endpoint
              end
              next if relation == :endpoint_touch && edge.polygon_index == other.polygon_index
              next if relation == :coincident &&
                      edge.polygon_index != other.polygon_index &&
                      edge.start_key == other.end_key

              return false
            end
          end
        end
        true
      end

      def segment_relation(first, second)
        a, b = first.start_key, first.end_key
        c, d = second.start_key, second.end_key
        cross1 = cross(a, b, c)
        cross2 = cross(a, b, d)
        cross3 = cross(c, d, a)
        cross4 = cross(c, d, b)

        if cross1.zero? && cross2.zero? && cross3.zero? && cross4.zero?
          return :coincident if geometry_key(a, b) == geometry_key(c, d)
          return :overlap if collinear_segments_overlap?(a, b, c, d)

          return :none
        end

        proper = ((cross1 * cross2) < 0) && ((cross3 * cross4) < 0)
        return :cross if proper

        endpoint_of_first = [a, b]
        endpoint_of_second = [c, d]
        return :cross if endpoint_of_first.any? { |key| point_on_segment?(key, point_from_key(c), point_from_key(d)) }
        return :cross if endpoint_of_second.any? { |key| point_on_segment?(key, point_from_key(a), point_from_key(b)) }

        return :endpoint_touch if (endpoint_of_first & endpoint_of_second).any?

        :none
      end

      def cross(start_key, end_key, point_key)
        (end_key[0] - start_key[0]) * (point_key[1] - start_key[1]) -
          (end_key[1] - start_key[1]) * (point_key[0] - start_key[0])
      end

      def collinear_segments_overlap?(a, b, c, d)
        axis = (a[0] - b[0]).abs > (a[1] - b[1]).abs ? 0 : 1
        left1, right1 = [a[axis], b[axis]].minmax
        left2, right2 = [c[axis], d[axis]].minmax
        [left1, left2].max < [right1, right2].min
      end

      def point_from_key(key)
        Point.from_x_y(key[0], key[1])
      end

      def apply_proposals(proposals)
        proposals.each do |ring_key, proposal|
          proposal[:retired].each { |id| @retired_edge_ids << id }
          proposal[:added].each do |edge|
            edge.id = next_edge_id
            @edges << edge
            @ring_edges[ring_key] << edge
          end
        end
      end

      def next_edge_id
        @next_edge_id += 1
      end

      def simplify_nonshared_rings
        @rings.each do |ring|
          key = ring_key(ring)
          edges = @ring_edges[key].reject { |edge| @retired_edge_ids.include?(edge.id) }
          cycle = cycle_keys(edges)
          next if cycle.empty?

          if ring.closed
            simplify_closed_nonshared_ring(ring, cycle)
          else
            simplify_open_nonshared_ring(ring, cycle)
          end
        end
      end

      def simplify_closed_nonshared_ring(ring, cycle)
        unique_cycle = cycle[0...-1]
        barriers = barrier_indices(unique_cycle)
        if barriers.length < 2
          centroid = centroid_key(unique_cycle)
          first = (0...unique_cycle.length).max_by do |index|
            distance_points(unique_cycle[index], centroid)
          end || 0
          second = (0...unique_cycle.length).reject { |index| index == first }.max_by do |index|
            distance_points(unique_cycle[index], unique_cycle[first])
          end || (first + 1) % unique_cycle.length
          barriers << first << second
        end
        barriers.sort!.uniq!
        barriers.each_with_index do |start_position, order|
          end_position = barriers[(order + 1) % barriers.length]
          next if start_position == end_position

          path = cyclic_path(unique_cycle, start_position, end_position)
          simplify_ring_path(ring, path)
        end
      end

      def simplify_open_nonshared_ring(ring, path)
        barriers = barrier_indices(path)
        barriers << 0 << path.length - 1
        barriers.sort!.uniq!
        barriers.each_cons(2) do |start_index, end_index|
          simplify_ring_path(ring, path[start_index..end_index])
        end
      end

      def barrier_indices(keys)
        (0...keys.length).select do |index|
          key = keys[index]
          @anchors.include?(key) || junction_node?(key) || shared_node?(key)
        end
      end

      def shared_node?(key)
        @shared_groups.any? do |group|
          group.start_key == key || group.end_key == key
        end
      end

      def cyclic_path(keys, start_position, end_position)
        path = []
          position = start_position
          loop do
            path << keys[position]
            break if position == end_position && path.length > 1

            position = (position + 1) % keys.length
            break if position == start_position && path.length > 1
          end
        path
      end

      def centroid_key(keys)
        [keys.sum { |key| key[0] } / keys.length.to_f,
         keys.sum { |key| key[1] } / keys.length.to_f]
      end

      def distance_points(first, second)
        Math.sqrt((first[0] - second[0])**2 + (first[1] - second[1])**2)
      end

      def simplify_ring_path(ring, path)
        return unless path.length > 2

        edges = directed_edges_for_path(ring, path)
        return unless edges.length == path.length - 1
        return if edges.any?(&:shared) || edges.any?(&:candidate)
        return if edges.any? { |edge| @retired_edge_ids.include?(edge.id) }

        simplify_nonshared_interval(ring, path, 0, path.length - 1, edges.map(&:id))
      end

      def directed_edges_for_path(ring, path)
        path.each_cons(2).map do |start_key, end_key|
          @ring_edges[ring_key(ring)].find do |edge|
            next false if @retired_edge_ids.include?(edge.id)

            edge.start_key == start_key && edge.end_key == end_key
          end
        end
      end

      def simplify_nonshared_interval(ring, path, start_index, end_index, edge_ids)
        return if end_index - start_index <= 1
        return unless (start_index...end_index).any? do |index|
          active_edge_between?(ring, path[index], path[index + 1])
        end
        if path[start_index] == path[end_index]
          retain_interval(path[start_index], path[end_index], 0)
          return
        end

        farthest_index = (start_index + 1...end_index).max_by do |index|
          distance_to_segment(path[index], path[start_index], path[end_index])
        end
        farthest_distance = distance_to_segment(path[farthest_index],
                                                path[start_index], path[end_index])

        if farthest_distance <= @epsilon
          proposals = nonshared_proposals(ring, path, start_index, end_index, edge_ids)
          if valid_proposals?(proposals)
            apply_proposals(proposals)
          else
            retain_interval(path[start_index], path[end_index], farthest_distance)
            simplify_nonshared_interval(ring, path, start_index, farthest_index, edge_ids)
            simplify_nonshared_interval(ring, path, farthest_index, end_index, edge_ids)
          end
        else
          simplify_nonshared_interval(ring, path, start_index, farthest_index, edge_ids)
          simplify_nonshared_interval(ring, path, farthest_index, end_index, edge_ids)
        end
      end

      def nonshared_proposals(ring, path, start_index, end_index, edge_ids)
        start_key = path[start_index]
        end_key = path[end_index]
        retired = path[start_index...end_index].each_cons(2).map do |segment_start, segment_end|
          edge_ids.find do |id|
            edge = @edges.find { |item| item.id == id && !@retired_edge_ids.include?(id) }
            edge&.start_key == segment_start && edge&.end_key == segment_end
          end
        end.compact
        edge = Edge.new(
          id: next_edge_id,
          polygon_index: ring.polygon_index,
          ring_index: ring.ring_index,
          start_key: start_key,
          end_key: end_key,
          start_point: @all_points[start_key],
          end_point: @all_points[end_key],
          shared: false,
          candidate: true
        )
        { ring_key(ring) => { retired: Set.new(retired), added: [edge] } }
      end

      def active_edge_between?(ring, start_key, end_key)
        @ring_edges[ring_key(ring)].any? do |edge|
          !@retired_edge_ids.include?(edge.id) &&
            edge.start_key == start_key && edge.end_key == end_key
        end
      end

      def build_result
        polygons = @polygons.map.with_index do |_polygon, polygon_index|
          rings = @polygons[polygon_index].rings.map.with_index do |_ring, ring_index|
            key = [polygon_index, ring_index]
            ring = @rings.find { |item| ring_key(item) == key }
            edges = @ring_edges[key].reject { |edge| @retired_edge_ids.include?(edge.id) }
            keys = cycle_keys(edges)
            keys << keys.first if ring.closed && keys.length > 1 && keys.first != keys.last
            points = keys.map { |point_key| @ring_point_maps[key][point_key] }
            LinearRing.from_points(points, ring.srid, ring.with_z, ring.with_m)
          end
          Polygon.from_linear_rings(rings, @polygons[polygon_index].srid,
                                   @polygons[polygon_index].with_z?,
                                   @polygons[polygon_index].with_m?)
        end
        multi_polygon = MultiPolygon.from_polygons(polygons)
        Result.new(multi_polygon: multi_polygon, warnings: @warnings)
      end
    end
  end
end
