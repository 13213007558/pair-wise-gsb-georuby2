require 'geo_ruby/simple_features/polygon'

module GeoRuby
  module SimpleFeatures
    # Jointly simplifies a set of polygons (for instance adjacent
    # administrative areas) so that no gaps or overlaps appear along
    # shared boundaries, which per-polygon LineString#simplify cannot
    # guarantee.
    #
    # Shared boundaries are detected as chains of edges whose coordinates
    # are exactly equal (in reverse order) in two polygons. The shared
    # edge may be split into a different number of collinear short edges
    # in each polygon: coincident chains are merged before simplification,
    # so callers do not have to harmonize the segmentation first.
    #
    # Guarantees:
    # * the same shared arc keeps the same vertices in every polygon;
    # * junction vertices (where the boundaries of two or more polygons
    #   meet) and user supplied anchor points are never removed;
    # * no new crossings or contacts are introduced, holes stay inside
    #   their polygon and originally disjoint faces stay disjoint;
    # * when a chain cannot be simplified within the tolerance without
    #   breaking those constraints, the original chain is kept and a
    #   reason is appended to #warnings.
    #
    # Usage:
    #   simplifier = GeoRuby::SimpleFeatures::JointSimplifier.new(polygons, 0.5)
    #   simplified = simplifier.simplify
    #   simplifier.warnings # => ["kept original boundary ...", ...]
    class JointSimplifier
      # Relative tolerance used to decide whether a point lies exactly
      # on a segment (shared boundaries are expected to use identical
      # coordinates, so this only absorbs floating point noise).
      ON_SEGMENT_TOLERANCE = 1e-9

      # One boundary chain between two junction vertices.
      Chain = Struct.new(:ring_index, :points, :final)

      # Reasons why some chains were left unsimplified (strings).
      attr_reader :warnings

      # +polygons+ is an array of Polygon (or a MultiPolygon), +epsilon+
      # the simplification tolerance (same semantics as
      # LineString#simplify) and +anchors+ an optional list of points
      # (Point instances or [x, y] pairs) that must be preserved.
      def initialize(polygons, epsilon = 1, anchors: [])
        @input_is_multi = polygons.is_a?(MultiPolygon)
        @polygons = (@input_is_multi ? polygons.geometries : polygons).to_a
        @epsilon = epsilon
        @anchors = anchors
        @warnings = []
      end

      # Simplifies every polygon jointly. Returns an array of Polygon
      # (a MultiPolygon when the input was one), in the input order.
      def simplify
        collect_rings
        detect_junctions
        process_anchors
        insert_junction_points
        split_into_chains
        group_chains
        simplify_groups
        rebuild_polygons
      end

      # Convenience class method when the warnings are not needed.
      def self.simplify(polygons, epsilon = 1, anchors: [])
        new(polygons, epsilon, anchors: anchors).simplify
      end

      private

      def coord_key(point)
        [point.x, point.y]
      end

      def cross(a, b, point)
        (b.x - a.x) * (point.y - a.y) - (b.y - a.y) * (point.x - a.x)
      end

      # Collects the rings of every polygon as unclosed point arrays.
      def collect_rings
        @rings = []
        @ring_polygon = []
        @ring_meta = []
        @polygons.each_with_index do |polygon, polygon_index|
          polygon.rings.each do |ring|
            points = ring.points.dup
            points.pop while points.size > 1 &&
                             coord_key(points.last) == coord_key(points.first)
            @rings << points
            @ring_polygon << polygon_index
            @ring_meta << [ring.srid, ring.with_z, ring.with_m]
          end
        end
      end

      # Junctions are vertices shared by at least two distinct polygons
      # where the boundary graph branches or ends: shared vertices in
      # the interior of a shared boundary have exactly two distinct
      # neighbours and are ordinary simplifiable points.
      def detect_junctions
        @junction_points = {}
        polygons_of = Hash.new { |hash, key| hash[key] = [] }
        neighbors_of = Hash.new { |hash, key| hash[key] = [] }
        @rings.each_with_index do |points, ring_index|
          points.each_with_index do |point, i|
            key = coord_key(point)
            owners = polygons_of[key]
            polygon_index = @ring_polygon[ring_index]
            owners << polygon_index unless owners.include?(polygon_index)
            neighbors = neighbors_of[key]
            previous = points[(i - 1) % points.size]
            following = points[(i + 1) % points.size]
            [previous, following].each do |other|
              other_key = coord_key(other)
              neighbors << other_key unless neighbors.include?(other_key)
            end
          end
        end
        polygons_of.each do |key, owners|
          next unless owners.size > 1 && neighbors_of[key].size != 2
          @junction_points[key] = find_vertex(key)
        end
      end

      def find_vertex(key)
        @rings.each do |points|
          found = points.find { |point| coord_key(point) == key }
          return found if found
        end
        nil
      end

      # Anchors matching an existing vertex become junctions; anchors
      # lying exactly on a boundary segment are inserted as new junction
      # vertices; anything else is reported and ignored.
      def process_anchors
        @anchors.each do |anchor|
          key = anchor.respond_to?(:x) ? [anchor.x, anchor.y] : [anchor[0], anchor[1]]
          vertex = find_vertex(key)
          if vertex
            @junction_points[key] ||= vertex
          elsif on_any_boundary?(key)
            @junction_points[key] ||= Point.from_x_y(key[0], key[1])
          else
            @warnings << "anchor (#{key[0]}, #{key[1]}) ignored: " \
                         'not located on any polygon boundary'
          end
        end
      end

      def on_any_boundary?(key)
        point = Point.from_x_y(key[0], key[1])
        @rings.any? do |points|
          points.each_with_index.any? do |a, i|
            on_segment?(a, points[(i + 1) % points.size], point)
          end
        end
      end

      # Inserts junction points that lie in the middle of an edge (for
      # instance a T-junction, or a shared edge split into a different
      # number of collinear edges in the neighbouring polygon).
      def insert_junction_points
        @rings.map! do |points|
          result = []
          points.each_with_index do |a, i|
            b = points[(i + 1) % points.size]
            result << a
            on_edge = @junction_points.values.select do |point|
              strictly_on_segment?(a, b, point)
            end
            on_edge.sort_by! { |point| segment_param(a, b, point) }
            result.concat(on_edge)
          end
          result
        end
      end

      def strictly_on_segment?(a, b, point)
        key = coord_key(point)
        key != coord_key(a) && key != coord_key(b) && on_segment?(a, b, point)
      end

      def on_segment?(a, b, point)
        key = coord_key(point)
        return key == coord_key(a) || key == coord_key(b) if coord_key(a) == coord_key(b)
        length = Math.hypot(b.x - a.x, b.y - a.y)
        tolerance = ON_SEGMENT_TOLERANCE * (length * length + 1.0)
        return false if cross(a, b, point).abs > tolerance
        margin = ON_SEGMENT_TOLERANCE * (length + 1.0)
        point.x >= [a.x, b.x].min - margin && point.x <= [a.x, b.x].max + margin &&
          point.y >= [a.y, b.y].min - margin && point.y <= [a.y, b.y].max + margin
      end

      def segment_param(a, b, point)
        dx = b.x - a.x
        dy = b.y - a.y
        ((point.x - a.x) * dx + (point.y - a.y) * dy) / (dx * dx + dy * dy)
      end

      # Splits every ring into chains between consecutive junctions.
      def split_into_chains
        @chains = []
        @ring_chains = Hash.new { |hash, key| hash[key] = [] }
        @rings.each_with_index do |points, ring_index|
          junctions = points.each_index.select do |i|
            @junction_points.key?(coord_key(points[i]))
          end
          if junctions.size < 2
            points = points.rotate(junctions.first) unless junctions.empty?
            add_closed_chains(ring_index, points)
          else
            junctions.each_with_index do |start, n|
              stop = junctions[(n + 1) % junctions.size]
              chain = []
              i = start
              loop do
                chain << points[i]
                break if i == stop
                i = (i + 1) % points.size
              end
              add_chain(ring_index, chain)
            end
          end
        end
      end

      # A ring with fewer than two junctions is split into two open
      # chains at the vertex farthest from its first point.
      def add_closed_chains(ring_index, points)
        farthest = (1...points.size).max_by do |i|
          points[0].euclidian_distance(points[i])
        end
        if farthest.nil?
          add_chain(ring_index, points + [points[0]])
        else
          add_chain(ring_index, points[0..farthest])
          add_chain(ring_index, points[farthest..-1] + [points[0]])
        end
      end

      def add_chain(ring_index, points)
        chain = Chain.new(ring_index, points, nil)
        @chains << chain
        @ring_chains[ring_index] << chain
      end

      # Groups chains connecting the same pair of junctions whose
      # vertices all lie on each other: coincident shared boundaries,
      # possibly segmented into a different number of collinear edges.
      def group_chains
        buckets = Hash.new { |hash, key| hash[key] = [] }
        @chains.each do |chain|
          key = [coord_key(chain.points.first), coord_key(chain.points.last)].sort
          buckets[key] << chain
        end
        @groups = []
        buckets.each_value do |bucket|
          bucket_groups = []
          bucket.each do |chain|
            group = bucket_groups.find { |members| coincident?(members.first, chain) }
            if group
              group << chain
            else
              bucket_groups << [chain]
            end
          end
          @groups.concat(bucket_groups)
        end
      end

      def coincident?(chain_a, chain_b)
        vertices_on_polyline?(chain_a.points, chain_b.points) &&
          vertices_on_polyline?(chain_b.points, chain_a.points)
      end

      def vertices_on_polyline?(vertices, polyline)
        segments = polyline.each_cons(2).to_a
        vertices.all? do |point|
          segments.any? { |a, b| on_segment?(a, b, point) }
        end
      end

      # Simplifies every group of coincident chains once, so that all
      # polygons sharing a boundary keep exactly the same vertices.
      def simplify_groups
        @groups.each do |members|
          merged = merge_group(members)
          simplified = douglas_peucker(merged)
          if simplified.size >= merged.size
            members.each { |chain| chain.final = chain.points }
          elsif valid_replacement?(members, simplified)
            members.each { |chain| chain.final = orient(simplified, chain) }
          else
            members.each { |chain| chain.final = chain.points }
          end
        end
      end

      # Merges coincident chains into a single polyline holding the
      # vertices of every member, so the tolerance also holds for
      # polygons whose shared edge was segmented differently.
      def merge_group(members)
        oriented = members.map { |chain| orient(chain.points, members.first) }
        merged = oriented.max_by(&:size).dup
        oriented.each do |points|
          points[1..-2].each do |point|
            next if merged.any? { |q| coord_key(q) == coord_key(point) }
            index = (0...merged.size - 1).find do |i|
              on_segment?(merged[i], merged[i + 1], point)
            end
            merged.insert(index + 1, point) if index
          end
        end
        merged
      end

      def orient(points, chain)
        if coord_key(chain.points.first) == coord_key(points.first)
          points
        else
          points.reverse
        end
      end

      # Ramer-Douglas-Peucker, same semantics as LineString#do_simplify.
      def douglas_peucker(points)
        return points.dup if points.size <= 2
        index = 0
        dmax = 0.0
        1.upto(points.size - 2) do |i|
          distance = orthogonal_distance(points[i], points.first, points.last)
          if distance > dmax
            index = i
            dmax = distance
          end
        end
        return [points.first, points.last] if dmax < @epsilon
        left = douglas_peucker(points[0..index])
        right = douglas_peucker(points[index..-1])
        left[0..-2] + right
      end

      # Orthogonal distance from +point+ to the segment ab, computed in
      # floating point (Point#orthogonal_distance suffers from integer
      # division when the coordinates are Integers).
      def orthogonal_distance(point, a, b)
        dx = b.x - a.x
        dy = b.y - a.y
        length_sq = dx * dx + dy * dy
        return point.euclidian_distance(a) if length_sq.zero?
        param = ((point.x - a.x) * dx + (point.y - a.y) * dy) / length_sq.to_f
        param = 0.0 if param.negative?
        param = 1.0 if param > 1.0
        Math.sqrt((point.x - a.x - param * dx)**2 + (point.y - a.y - param * dy)**2)
      end

      # A simplified chain is accepted only if it does not intersect
      # itself, keeps every affected ring non-degenerate and does not
      # cross or touch any other chain (except at shared junction
      # endpoints). Otherwise the original chain is kept and the reason
      # is recorded.
      def valid_replacement?(members, simplified)
        description = describe(members)
        if self_intersecting?(simplified)
          @warnings << "kept original boundary #{description}: " \
                       'simplification would self-intersect'
          return false
        end
        members.map(&:ring_index).uniq.each do |ring_index|
          next if ring_valid_after?(ring_index, members, simplified)
          @warnings << "kept original boundary #{description}: simplification " \
                       "would degenerate a ring of polygon #{@ring_polygon[ring_index]}"
          return false
        end
        @chains.each do |other|
          next if members.include?(other)
          contact = forbidden_contact(simplified, other)
          next unless contact
          polygon_index = @ring_polygon[other.ring_index]
          @warnings << "kept original boundary #{description}: simplification " \
                       "would cross or touch polygon #{polygon_index} " \
                       "near (#{contact[0]}, #{contact[1]})"
          return false
        end
        true
      end

      def describe(members)
        first = members.first.points.first
        last = members.first.points.last
        "(#{first.x}, #{first.y}) -> (#{last.x}, #{last.y})"
      end

      def self_intersecting?(points)
        segments = points.each_cons(2).to_a
        segments.each_with_index do |seg_a, i|
          segments.each_with_index do |seg_b, j|
            next if j <= i + 1
            return true if segment_contact(seg_a[0], seg_a[1], seg_b[0], seg_b[1])
          end
        end
        false
      end

      # Returns the coordinate of the first contact between the
      # simplified chain and another chain that is not a shared
      # junction endpoint, or nil when there is none.
      def forbidden_contact(simplified, other)
        other_points = other.final || other.points
        allowed = [coord_key(other_points.first), coord_key(other_points.last)] &
                  [coord_key(simplified.first), coord_key(simplified.last)]
        simplified.each_cons(2) do |a, b|
          other_points.each_cons(2) do |c, d|
            contact = segment_contact(a, b, c, d)
            return contact if contact && !allowed.include?(contact)
          end
        end
        nil
      end

      # Returns the [x, y] contact of two segments (any crossing,
      # touch or collinear overlap) or nil when they are disjoint.
      def segment_contact(a, b, c, d)
        return point_segment_contact(a, c, d) if coord_key(a) == coord_key(b)
        return point_segment_contact(c, a, b) if coord_key(c) == coord_key(d)
        o1 = cross(a, b, c)
        o2 = cross(a, b, d)
        o3 = cross(c, d, a)
        o4 = cross(c, d, b)
        if o1.zero? && o2.zero? && o3.zero? && o4.zero?
          collinear_contact(a, b, c, d)
        elsif (o1 * o2).negative? && (o3 * o4).negative?
          line_intersection(a, b, c, d)
        else
          return coord_key(c) if o1.zero? && on_segment?(a, b, c)
          return coord_key(d) if o2.zero? && on_segment?(a, b, d)
          return coord_key(a) if o3.zero? && on_segment?(c, d, a)
          return coord_key(b) if o4.zero? && on_segment?(c, d, b)
          nil
        end
      end

      def point_segment_contact(point, c, d)
        on_segment?(c, d, point) ? coord_key(point) : nil
      end

      def collinear_contact(a, b, c, d)
        if (b.x - a.x).abs >= (b.y - a.y).abs
          lo = [[a.x, b.x].min, [c.x, d.x].min].max
          hi = [[a.x, b.x].max, [c.x, d.x].max].min
          return nil if lo > hi
          # a positive-length overlap is reported through its midpoint,
          # which can never be confused with a shared junction endpoint
          x = hi == lo ? lo : (lo + hi) / 2.0
          [x, a.y + (x - a.x) * (b.y - a.y) / (b.x - a.x)]
        else
          lo = [[a.y, b.y].min, [c.y, d.y].min].max
          hi = [[a.y, b.y].max, [c.y, d.y].max].min
          return nil if lo > hi
          y = hi == lo ? lo : (lo + hi) / 2.0
          [a.x + (y - a.y) * (b.x - a.x) / (b.y - a.y), y]
        end
      end

      def line_intersection(a, b, c, d)
        t = ((c.x - a.x) * (d.y - c.y) - (c.y - a.y) * (d.x - c.x)) /
            ((b.x - a.x) * (d.y - c.y) - (b.y - a.y) * (d.x - c.x))
        [a.x + t * (b.x - a.x), a.y + t * (b.y - a.y)]
      end

      # Simulates a ring after replacing its chains in the group and
      # checks it still has at least 3 distinct vertices and a
      # non-zero area.
      def ring_valid_after?(ring_index, members, simplified)
        replaced = {}
        members.each do |chain|
          replaced[chain] = orient(simplified, chain) if chain.ring_index == ring_index
        end
        points = ring_points_after(ring_index, replaced)
        return false if points.map { |point| coord_key(point) }.uniq.size < 3
        scale = points.map { |point| [point.x.abs, point.y.abs].max }.max || 0.0
        shoelace(points).abs > ON_SEGMENT_TOLERANCE * (scale * scale + 1.0)
      end

      def ring_points_after(ring_index, replaced)
        points = []
        @ring_chains[ring_index].each do |chain|
          sequence = replaced[chain] || chain.final || chain.points
          points.concat(sequence[0..-2])
        end
        points
      end

      def shoelace(points)
        sum = 0.0
        points.each_with_index do |a, i|
          b = points[(i + 1) % points.size]
          sum += a.x * b.y - b.x * a.y
        end
        sum / 2.0
      end

      def rebuild_polygons
        rings_by_polygon = Hash.new { |hash, key| hash[key] = [] }
        @ring_polygon.each_with_index do |polygon_index, ring_index|
          rings_by_polygon[polygon_index] << ring_index
        end
        polygons = @polygons.each_with_index.map do |polygon, polygon_index|
          linear_rings = rings_by_polygon[polygon_index].map do |ring_index|
            points = []
            @ring_chains[ring_index].each { |chain| points.concat(chain.final[0..-2]) }
            points << points.first
            srid, with_z, with_m = @ring_meta[ring_index]
            LinearRing.from_points(points, srid, with_z, with_m)
          end
          Polygon.from_linear_rings(linear_rings, polygon.srid,
                                    polygon.with_z, polygon.with_m)
        end
        @input_is_multi ? MultiPolygon.from_polygons(polygons) : polygons
      end
    end
  end
end
