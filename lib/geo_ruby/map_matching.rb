module GeoRuby
  # Offline map matching for time-ordered GPS observations against a directed
  # road polyline network.
  #
  # The matcher is an HMM solved jointly with Viterbi (Newson & Krumm style):
  #
  # * Emission cost of a candidate is gaussian: d^2 / (2 * sigma^2), where d is
  #   the perpendicular distance between the observation and its projection on
  #   the edge. Candidates are projections anywhere on a segment, not only on
  #   road vertices.
  # * Transition cost compares the network route distance between two
  #   consecutive candidates with the straight-line displacement of the two
  #   observations: |route_length - displacement| / beta. The route respects
  #   the directed graph, so unreachable transitions cost infinity.
  #
  # Defined semantics:
  #
  # * One-way edges can only be traversed from their first coordinate to their
  #   last. Moving backwards on a one-way edge requires a loop through the
  #   rest of the network; if no such loop exists the transition is impossible.
  # * On a two-way edge, moving backwards along the same edge is allowed and
  #   costs the absolute offset difference.
  # * Disconnected regions have no route between them; the transition is
  #   impossible. When no candidate of an observation is reachable from the
  #   previous one, the sequence is split into explicit segments. Segments are
  #   never joined by a synthetic straight line.
  # * Observations with a nil point (missing fixes) or with no candidate inside
  #   +max_radius+ are skipped and reported in Result#unmatched. The transition
  #   that bridges a gap uses the displacement between the surrounding valid
  #   observations.
  # * Ties are broken deterministically: lower total cost, then lower
  #   perpendicular distance, then smaller edge id (compared as strings), then
  #   smaller offset along the edge.
  #
  # Coordinates are planar (e.g. meters in a projected CRS); all distances and
  # parameters share the same unit.
  module MapMatching
    EPSILON = 1e-9

    # Normalizes a point given as [x, y] or as an object answering x/y.
    def self.normalize_point(point)
      if point.respond_to?(:x) && point.respond_to?(:y)
        [point.x.to_f, point.y.to_f]
      else
        [point[0].to_f, point[1].to_f]
      end
    end

    # A directed (or two-way) road polyline.
    #
    # +nodes+ optionally gives explicit endpoint node ids. By default endpoints
    # are identified by their exact coordinates, so crossing polylines connect
    # only when they share a vertex. Use explicit node ids to model overpasses
    # that cross geometrically without connecting.
    class Edge
      attr_reader :id, :coordinates, :oneway
      attr_accessor :from_node, :to_node

      def initialize(id, coordinates, oneway: false, nodes: nil)
        @id = id
        @coordinates = coordinates.map { |c| MapMatching.normalize_point(c) }
        raise ArgumentError, 'an edge needs at least 2 coordinates' if @coordinates.size < 2

        @oneway = oneway
        @node_ids = nodes
        @cum = [0.0]
        (0...@coordinates.size - 1).each do |i|
          ax, ay = @coordinates[i]
          bx, by = @coordinates[i + 1]
          @cum << @cum[-1] + Math.hypot(bx - ax, by - ay)
        end
      end

      def length
        @cum[-1]
      end

      def from_node_key
        @node_ids ? @node_ids[0] : @coordinates.first
      end

      def to_node_key
        @node_ids ? @node_ids[1] : @coordinates.last
      end

      # Projects a point on the polyline. Returns the offset along the edge,
      # the projected point and the perpendicular distance.
      def project(point)
        px, py = point
        best = nil
        (0...@coordinates.size - 1).each do |i|
          ax, ay = @coordinates[i]
          bx, by = @coordinates[i + 1]
          dx = bx - ax
          dy = by - ay
          seg_len_sq = (dx * dx) + (dy * dy)
          t = seg_len_sq.zero? ? 0.0 : (((px - ax) * dx) + ((py - ay) * dy)) / seg_len_sq
          t = 0.0 if t.negative?
          t = 1.0 if t > 1.0
          qx = ax + (t * dx)
          qy = ay + (t * dy)
          distance = Math.hypot(px - qx, py - qy)
          offset = @cum[i] + (t * Math.sqrt(seg_len_sq))
          if best.nil? || distance < best[:distance] - EPSILON ||
             ((distance - best[:distance]).abs <= EPSILON && offset < best[:offset])
            best = { offset: offset, point: [qx, qy], distance: distance }
          end
        end
        best
      end

      # Interpolates the point at the given offset along the edge.
      def point_at(offset)
        offset = 0.0 if offset.negative?
        offset = length if offset > length
        (0...@coordinates.size - 1).each do |i|
          seg_len = @cum[i + 1] - @cum[i]
          next if @cum[i + 1] < offset && i < @coordinates.size - 2

          t = seg_len.zero? ? 0.0 : (offset - @cum[i]) / seg_len
          ax, ay = @coordinates[i]
          bx, by = @coordinates[i + 1]
          return [ax + (t * (bx - ax)), ay + (t * (by - ay))]
        end
        @coordinates.last
      end
    end

    # A candidate matched position: the projection of an observation on an edge.
    class Candidate
      attr_reader :edge, :offset, :point, :distance

      def initialize(edge, offset, point, distance)
        @edge = edge
        @offset = offset
        @point = point
        @distance = distance
      end

      def edge_id
        @edge.id
      end
    end

    # One traversal of one edge inside a Route. +from_offset+ > +to_offset+
    # means the edge is traversed against its digitized direction (only
    # possible on two-way edges).
    Leg = Struct.new(:edge_id, :from_offset, :to_offset)

    # A network route between two candidates.
    class Route
      attr_reader :length, :legs

      def initialize(length, legs)
        @length = length
        @legs = legs
      end

      def edge_ids
        @legs.map(&:edge_id)
      end
    end

    Arc = Struct.new(:to, :weight, :edge_id, :from_offset, :to_offset)

    # A directed road network built from Edge objects.
    class Network
      attr_reader :edges

      def initialize(edges)
        @edges = edges
        @by_id = {}
        @arcs = Hash.new { |hash, key| hash[key] = [] }
        edges.each do |edge|
          raise ArgumentError, "duplicate edge id: #{edge.id}" if @by_id.key?(edge.id)

          @by_id[edge.id] = edge
          edge.from_node = edge.from_node_key
          edge.to_node = edge.to_node_key
          @arcs[edge.from_node] << Arc.new(edge.to_node, edge.length, edge.id, 0.0, edge.length)
          unless edge.oneway
            @arcs[edge.to_node] << Arc.new(edge.from_node, edge.length, edge.id, edge.length, 0.0)
          end
        end
      end

      # All projections of +point+ within +max_radius+, sorted deterministically.
      def candidates(point, max_radius)
        p = MapMatching.normalize_point(point)
        found = []
        @edges.each do |edge|
          proj = edge.project(p)
          next unless proj[:distance] <= max_radius

          found << Candidate.new(edge, proj[:offset], proj[:point], proj[:distance])
        end
        found.sort_by { |c| [c.distance, c.edge.id.to_s, c.offset] }
      end

      # Shortest legal route between two candidates, or nil when the network
      # cannot connect them (disconnected regions, one-way traps).
      def route(from, to)
        options = []
        if from.edge.equal?(to.edge)
          # Moving forward along the edge always follows the allowed direction.
          if to.offset >= from.offset - EPSILON
            options << Route.new((to.offset - from.offset).abs,
                                 [Leg.new(from.edge.id, from.offset, to.offset)])
          end
          # Moving backwards along the same edge is only allowed on two-way edges.
          unless from.edge.oneway
            options << Route.new((to.offset - from.offset).abs,
                                 [Leg.new(from.edge.id, from.offset, to.offset)])
          end
        end
        exits = [[from.edge.length - from.offset, from.edge.to_node, from.edge.length]]
        exits << [from.offset, from.edge.from_node, 0.0] unless from.edge.oneway
        entries = [[to.offset, to.edge.from_node, 0.0]]
        entries << [to.edge.length - to.offset, to.edge.to_node, to.edge.length] unless to.edge.oneway
        exits.each do |exit_cost, exit_node, exit_offset|
          entries.each do |entry_cost, entry_node, entry_offset|
            path = shortest_path(exit_node, entry_node)
            next unless path

            middle_cost, middle_arcs = path
            legs = [Leg.new(from.edge.id, from.offset, exit_offset)]
            middle_arcs.each { |arc| legs << Leg.new(arc.edge_id, arc.from_offset, arc.to_offset) }
            legs << Leg.new(to.edge.id, entry_offset, to.offset)
            options << Route.new(exit_cost + middle_cost + entry_cost, legs)
          end
        end
        return nil if options.empty?

        options.min_by { |r| [r.length, r.legs.map { |l| l.edge_id.to_s }.join('/')] }
      end

      # Dijkstra shortest path between two node keys.
      # Returns [cost, [Arc, ...]] or nil when unreachable.
      def shortest_path(from_node, to_node)
        return [0.0, []] if from_node == to_node

        dist = { from_node => 0.0 }
        prev = {}
        visited = {}
        queue = [[0.0, from_node]]
        until queue.empty?
          index = queue.each_index.min_by { |i| queue[i][0] }
          cost, node = queue.delete_at(index)
          next if visited[node]

          visited[node] = true
          return [cost, unwind(prev, from_node, to_node)] if node == to_node

          @arcs[node].each do |arc|
            new_cost = cost + arc.weight
            next unless dist[arc.to].nil? || new_cost < dist[arc.to] - EPSILON

            dist[arc.to] = new_cost
            prev[arc.to] = [node, arc]
            queue << [new_cost, arc.to]
          end
        end
        nil
      end

      private

      def unwind(prev, from_node, to_node)
        arcs = []
        node = to_node
        while node != from_node
          parent, arc = prev[node]
          arcs.unshift(arc)
          node = parent
        end
        arcs
      end
    end

    # A matched position on the road network.
    MatchedPoint = Struct.new(:observation_index, :time, :edge_id, :offset,
                              :point, :observation, :distance, :emission_cost)

    # A maximal run of consecutive observations matched without a network break.
    class Segment
      attr_reader :matched_points, :routes, :transition_costs

      def initialize
        @matched_points = []
        @routes = []
        @transition_costs = []
      end

      def emission_cost
        @matched_points.sum(&:emission_cost)
      end

      def transition_cost
        @transition_costs.sum
      end

      # Total cost of the segment: emissions plus transitions.
      def score
        emission_cost + transition_cost
      end

      def observation_range
        return nil if @matched_points.empty?

        @matched_points.first.observation_index..@matched_points.last.observation_index
      end
    end

    # The result of a Matcher#match call.
    class Result
      attr_reader :segments, :unmatched

      def initialize(segments, unmatched)
        @segments = segments
        @unmatched = unmatched
      end

      def matched_points
        @segments.flat_map(&:matched_points)
      end
    end

    # Joint (Viterbi) matcher over a whole observation sequence.
    class Matcher
      attr_reader :network, :sigma, :beta, :max_radius

      def initialize(network, sigma: 10.0, beta: 50.0, max_radius: 50.0)
        @network = network
        @sigma = sigma
        @beta = beta
        @max_radius = max_radius
        @route_cache = {}
      end

      # +observations+ is a time-ordered array of Hashes/structs with :time and
      # :point, or [time, point] pairs. A nil point marks a missing observation.
      def match(observations)
        lattice = []
        unmatched = []
        observations.each_with_index do |obs, index|
          time, point = normalize_observation(obs)
          if point.nil?
            unmatched << index
            next
          end
          found = @network.candidates(point, @max_radius)
          if found.empty?
            unmatched << index
          else
            lattice << { index: index, time: time, point: point, candidates: found }
          end
        end
        return Result.new([], unmatched) if lattice.empty?

        segments = blocks(lattice).map { |block| viterbi(block) }
        Result.new(segments, unmatched)
      end

      private

      def normalize_observation(obs)
        time, point =
          if obs.is_a?(Hash)
            [obs[:time], obs[:point]]
          elsif obs.respond_to?(:time) && obs.respond_to?(:point)
            [obs.time, obs.point]
          else
            [obs[0], obs[1]]
          end
        [time, point && MapMatching.normalize_point(point)]
      end

      def emission_cost(candidate)
        (candidate.distance**2) / (2.0 * (@sigma**2))
      end

      def route_between(from, to)
        key = [from.object_id, to.object_id]
        @route_cache[key] = @network.route(from, to) unless @route_cache.key?(key)
        @route_cache[key]
      end

      # Splits the lattice into maximal blocks where every consecutive pair of
      # steps keeps at least one reachable candidate transition. A boundary is
      # an explicit segment break: blocks are never connected afterwards.
      def blocks(lattice)
        live = lattice.first[:candidates]
        result = []
        current = [lattice.first]
        (1...lattice.size).each do |k|
          step = lattice[k]
          reachable = step[:candidates].select do |cand|
            live.any? { |prev| route_between(prev, cand) }
          end
          if reachable.empty?
            result << current
            current = [step]
            live = step[:candidates]
          else
            current << step.merge(candidates: reachable)
            live = reachable
          end
        end
        result << current
        result
      end

      def viterbi(block)
        states = [block.first[:candidates].map do |cand|
          { cost: emission_cost(cand), prev: nil, route: nil, transition_cost: 0.0 }
        end]
        (1...block.size).each do |k|
          displacement = Math.hypot(block[k][:point][0] - block[k - 1][:point][0],
                                    block[k][:point][1] - block[k - 1][:point][1])
          states << block[k][:candidates].map do |cand|
            best = nil
            block[k - 1][:candidates].each_with_index do |prev_cand, pi|
              route = route_between(prev_cand, cand)
              next unless route

              transition = (route.length - displacement).abs / @beta
              arrival = states[k - 1][pi][:cost] + transition
              key = [arrival, prev_cand.distance, prev_cand.edge.id.to_s, prev_cand.offset]
              if best.nil? || (key <=> best[0]).negative?
                best = [key,
                        { cost: arrival + emission_cost(cand), prev: pi,
                          route: route, transition_cost: transition }]
              end
            end
            # Reachable candidates always have at least one finite predecessor.
            best[1]
          end
        end

        last = block.size - 1
        best_index = block[last][:candidates].each_with_index.min_by do |cand, ci|
          [states[last][ci][:cost], cand.distance, cand.edge.id.to_s, cand.offset]
        end[1]

        segment = Segment.new
        ci = best_index
        (block.size - 1).downto(0) do |k|
          cand = block[k][:candidates][ci]
          state = states[k][ci]
          segment.matched_points.unshift(MatchedPoint.new(block[k][:index], block[k][:time],
                                                          cand.edge_id, cand.offset, cand.point,
                                                          block[k][:point], cand.distance,
                                                          emission_cost(cand)))
          if state[:route]
            segment.routes.unshift(state[:route])
            segment.transition_costs.unshift(state[:transition_cost])
          end
          ci = state[:prev] if state[:prev]
        end
        segment
      end
    end
  end
end

