# frozen_string_literal: true

module GeoRuby
  # Offline map matching: snaps a time-ordered sequence of noisy (and
  # possibly missing) GPS observations onto a directed road network and
  # recovers the route driven between consecutive matched positions.
  #
  # The matcher is a Hidden Markov Model in the style of Newson and Krumm
  # (2009), decoded jointly over the whole sequence with Viterbi:
  #
  # * Emission: Gaussian on the distance between an observation and its
  #   projection on a road edge (+sigma+).
  # * Transition: exponential on the absolute difference between the
  #   network route distance (reachable and direction-aware) and the
  #   straight-line displacement of the two observations (+beta+).
  #
  # Defined behaviours:
  #
  # * Candidates are projections anywhere on an edge polyline, not only
  #   on road vertices.
  # * One-way edges are only traversed in their polyline direction.
  # * Moving backwards along a two-way edge is allowed and costs the
  #   travelled distance plus +u_turn_penalty+; on a one-way edge it is
  #   only possible by looping back through the network (impossible when
  #   no such loop exists).
  # * Edges connect only at shared vertices (identical coordinates).
  #   Crossing polylines without a shared vertex (e.g. an overpass) and
  #   disconnected components simply have no route between them.
  # * When no reachable transition exists between two consecutive
  #   observations the sequence is split into explicit segments (see
  #   MatchResult#breaks); segments are never joined by fake straight
  #   lines.
  # * Ties between equally scored paths are broken deterministically:
  #   candidates are ordered by [distance, edge id, fraction] and, among
  #   scores equal within 1e-9, the earliest candidate wins.
  #
  # Every segment reports its total log-probability score plus the
  # emission and transition parts.
  #
  # Coordinates must be planar (use a projected CRS, units in meters).
  module MapMatcher
    # Tolerance used when comparing scores for deterministic tie-breaking.
    SCORE_EPS = 1e-9
    SQRT_2PI = Math.sqrt(2.0 * Math::PI)

    # A single GPS observation. +x+/+y+ may be nil to represent a missing
    # fix; missing observations are skipped by the matcher (the chain
    # bridges over them) and reported as unmatched. +t+ is informational.
    Observation = Struct.new(:t, :x, :y) do
      def missing?
        x.nil? || y.nil?
      end
    end

    # A directed road edge. +points+ is an array of [x, y] pairs ordered
    # along the driving direction.
    class Edge
      attr_reader :id, :points, :length

      def initialize(id, points, oneway)
        @id = id
        @points = points.map { |p| [Float(p[0]), Float(p[1])] }
        @oneway = oneway ? true : false
        @cum = [0.0]
        (0...@points.length - 1).each do |i|
          @cum << @cum[-1] + segment_length(i)
        end
        @length = @cum[-1]
        raise ArgumentError, "edge #{id.inspect} has zero length" if @length <= 0
      end

      def oneway?
        @oneway
      end

      # Projects (px, py) on the polyline.
      # Returns [fraction (0..1), distance, [qx, qy]].
      def project(px, py)
        best = nil
        (0...@points.length - 1).each do |i|
          ax, ay = @points[i]
          bx, by = @points[i + 1]
          dx = bx - ax
          dy = by - ay
          len2 = dx * dx + dy * dy
          t = len2.zero? ? 0.0 : ((px - ax) * dx + (py - ay) * dy) / len2
          t = 0.0 if t.negative?
          t = 1.0 if t > 1.0
          qx = ax + t * dx
          qy = ay + t * dy
          d = Math.hypot(px - qx, py - qy)
          next if best && d >= best[0]

          along = @cum[i] + Math.sqrt(len2) * t
          best = [d, along / @length, [qx, qy]]
        end
        [best[1], best[0], best[2]]
      end

      # Point at +fraction+ along the polyline.
      def point_at(fraction)
        target = fraction * @length
        i = 0
        i += 1 while i < @cum.length - 2 && @cum[i + 1] < target
        seg = @cum[i + 1] - @cum[i]
        t = seg.zero? ? 0.0 : (target - @cum[i]) / seg
        ax, ay = @points[i]
        bx, by = @points[i + 1]
        [ax + t * (bx - ax), ay + t * (by - ay)]
      end

      private

      def segment_length(i)
        Math.hypot(@points[i + 1][0] - @points[i][0],
                   @points[i + 1][1] - @points[i][1])
      end
    end
    # A candidate match: the projection of an observation on an edge.
    class Candidate
      attr_accessor :emission_logprob
      attr_reader :edge, :fraction, :distance, :point

      def initialize(edge, fraction, distance, point)
        @edge = edge
        @fraction = fraction
        @distance = distance
        @point = point
      end
    end

    # A route along the network between two positions on edges.
    class Route
      # One polyline traversal: edge_id from from_fraction to to_fraction.
      # to_fraction < from_fraction means driving the edge backwards
      # (only possible on two-way edges).
      Step = Struct.new(:edge_id, :from_fraction, :to_fraction)

      attr_reader :steps, :length

      def initialize(steps, length)
        @steps = steps
        @length = length
      end

      def edge_ids
        @steps.map(&:edge_id).uniq
      end
    end

    # A directed network of road edges built from hashes like
    # { id:, points: [[x, y], ...], oneway: bool } or Edge instances.
    class RoadNetwork
      attr_reader :edges

      def initialize(edges)
        @edges = edges.map do |e|
          e.is_a?(Edge) ? e : Edge.new(e[:id], e[:points], e[:oneway])
        end
        @vertex_ids = {}
        @adjacency = Hash.new { |h, k| h[k] = [] }
        @start_vid = {}
        @end_vid = {}
        @edges.each do |edge|
          sv = vertex_id(edge.points.first)
          ev = vertex_id(edge.points.last)
          @start_vid[edge] = sv
          @end_vid[edge] = ev
          @adjacency[sv] << [ev, edge, true]
          @adjacency[ev] << [sv, edge, false] unless edge.oneway?
        end
      end

      # Candidate projections of (x, y) on every edge within +radius+,
      # ordered by [distance, edge id, fraction], at most +max_candidates+.
      def candidates(x, y, radius, max_candidates)
        found = []
        @edges.each do |edge|
          fraction, distance, point = edge.project(x, y)
          next unless distance <= radius

          found << Candidate.new(edge, fraction, distance, point)
        end
        found.sort_by! { |c| [c.distance, c.edge.id.to_s, c.fraction] }
        found.first(max_candidates)
      end

      # Shortest direction-aware route between two edge positions, or nil
      # when unreachable. +from+/+to+ respond to #edge and #fraction.
      def route(from, to, u_turn_penalty: 0.0)
        options = []
        if from.edge.equal?(to.edge)
          direct = same_edge_route(from, to, u_turn_penalty)
          options << direct if direct
        end
        net = network_route(from, to)
        options << net if net
        options.min_by(&:length)
      end

      private

      def vertex_id(point)
        key = [(point[0] * 1e6).round, (point[1] * 1e6).round].join(":")
        @vertex_ids.fetch(key) { @vertex_ids[key] = @vertex_ids.length }
      end

      def same_edge_route(from, to, u_turn_penalty)
        edge = from.edge
        df = to.fraction - from.fraction
        if df >= 0
          Route.new([Route::Step.new(edge.id, from.fraction, to.fraction)],
                    df * edge.length)
        elsif !edge.oneway?
          Route.new([Route::Step.new(edge.id, from.fraction, to.fraction)],
                    -df * edge.length + u_turn_penalty)
        end
      end

      def network_route(from, to)
        dist = {}
        prev = {}
        add_start(dist, @end_vid[from.edge],
                  (1.0 - from.fraction) * from.edge.length)
        unless from.edge.oneway?
          add_start(dist, @start_vid[from.edge],
                    from.fraction * from.edge.length)
        end
        dijkstra(dist, prev)
        best = best_target(from, to, dist, prev)
        best && build_route(from, to, best, prev)
      end

      def add_start(dist, vid, cost)
        dist[vid] = cost if dist[vid].nil? || cost < dist[vid]
      end

      def dijkstra(dist, prev)
        settled = {}
        loop do
          vid = nil
          best = nil
          dist.each do |v, c|
            next if settled[v]
            next if best && c >= best

            vid = v
            best = c
          end
          break unless vid

          settled[vid] = true
          @adjacency[vid].each do |to_vid, edge, forward|
            nc = best + edge.length
            next if settled[to_vid]
            next if dist[to_vid] && nc >= dist[to_vid] - SCORE_EPS

            dist[to_vid] = nc
            prev[to_vid] = [vid, edge, forward]
          end
        end
      end

      def best_target(from, to, dist, prev)
        targets = [[@start_vid[to.edge], to.fraction * to.edge.length, true]]
        unless to.edge.oneway?
          targets << [@end_vid[to.edge],
                      (1.0 - to.fraction) * to.edge.length, false]
        end
        best = nil
        targets.each do |vid, extra, forward|
          next unless dist[vid]
          # Skip degenerate exit-and-reenter paths on the same edge.
          next if from.edge.equal?(to.edge) && !prev.key?(vid)

          total = dist[vid] + extra
          best = [total, vid, forward] if best.nil? || total < best[0] - SCORE_EPS
        end
        best
      end

      def build_route(from, to, best, prev)
        total, vid, enter_forward = best
        arcs = []
        while (p = prev[vid])
          arcs.unshift(p)
          vid = p[0]
        end
        exit_forward = vid == @end_vid[from.edge]
        steps = [Route::Step.new(from.edge.id, from.fraction,
                                 exit_forward ? 1.0 : 0.0)]
        arcs.each do |_from_vid, edge, forward|
          steps << Route::Step.new(edge.id, forward ? 0.0 : 1.0,
                                   forward ? 1.0 : 0.0)
        end
        steps << Route::Step.new(to.edge.id, enter_forward ? 0.0 : 1.0,
                                 to.fraction)
        Route.new(steps, total)
      end
    end
    Unmatched = Struct.new(:index, :reason)
    Match = Struct.new(:obs_index, :edge_id, :fraction, :point, :distance,
                       :emission_logprob)
    Transition = Struct.new(:from_index, :to_index, :route, :logprob)
    Segment = Struct.new(:start_index, :end_index, :matches, :transitions,
                         :score, :emission_score, :transition_score)
    MatchResult = Struct.new(:segments, :unmatched, :breaks)
    Entry = Struct.new(:candidate, :score, :back, :route, :trans_logprob,
                       :obs_index, keyword_init: true)

    # Joint (Viterbi) matcher. See the module docs for the model.
    class Matcher
      attr_reader :network, :sigma, :beta, :search_radius, :max_candidates,
                  :u_turn_penalty

      def initialize(network, sigma:, beta:, search_radius: nil,
                     max_candidates: 8, u_turn_penalty: 0.0)
        @network = network
        @sigma = Float(sigma)
        @beta = Float(beta)
        @search_radius = search_radius ? Float(search_radius) : 4.0 * @sigma
        @max_candidates = max_candidates
        @u_turn_penalty = Float(u_turn_penalty)
        raise ArgumentError, "sigma must be positive" unless @sigma.positive?
        raise ArgumentError, "beta must be positive" unless @beta.positive?
      end

      def emission_logprob(distance)
        d = distance / @sigma
        -0.5 * d * d - Math.log(@sigma * SQRT_2PI)
      end

      def transition_logprob(route_length, linear_distance)
        -(route_length - linear_distance).abs / @beta - Math.log(@beta)
      end

      # Accepts Observation, { t:, x:, y: } hashes or [x, y] / [t, x, y]
      # arrays. Returns a MatchResult.
      def match(observations)
        segments = []
        unmatched = []
        breaks = []
        layer = nil
        last_obs = nil
        seg_start = nil
        close = lambda do
          segments << build_segment(seg_start, layer) if layer
          layer = nil
          last_obs = nil
        end
        observations.each_with_index do |raw, index|
          obs = Matcher.normalize(raw)
          if obs.missing?
            unmatched << Unmatched.new(index, :missing)
            next
          end
          cands = @network.candidates(obs.x, obs.y, @search_radius,
                                      @max_candidates)
          cands.each { |c| c.emission_logprob = emission_logprob(c.distance) }
          if cands.empty?
            close.call
            unmatched << Unmatched.new(index, :no_candidate)
            next
          end
          if layer.nil?
            seg_start = index
            layer = seed_layer(cands, index)
          else
            linear = Math.hypot(obs.x - last_obs.x, obs.y - last_obs.y)
            new_layer = best_entries(layer, cands, linear, index)
            if new_layer.empty?
              breaks << index
              close.call
              seg_start = index
              layer = seed_layer(cands, index)
            else
              layer = new_layer
            end
          end
          last_obs = obs
        end
        close.call
        MatchResult.new(segments, unmatched, breaks)
      end

      def self.normalize(raw)
        case raw
        when Observation then raw
        when Hash then Observation.new(raw[:t], raw[:x], raw[:y])
        when Array
          raw.length == 2 ? Observation.new(nil, raw[0], raw[1])
                          : Observation.new(raw[0], raw[1], raw[2])
        else
          raise ArgumentError, "unsupported observation: " + raw.inspect
        end
      end

      private

      def seed_layer(cands, index)
        cands.map do |c|
          Entry.new(candidate: c, score: c.emission_logprob, obs_index: index)
        end
      end

      def best_entries(layer, cands, linear, obs_index)
        entries = []
        cands.each do |cand|
          best = nil
          layer.each do |entry|
            route = @network.route(entry.candidate, cand,
                                   u_turn_penalty: @u_turn_penalty)
            next unless route

            tlog = transition_logprob(route.length, linear)
            score = entry.score + tlog + cand.emission_logprob
            next if best && score <= best.score + SCORE_EPS

            best = Entry.new(candidate: cand, score: score, back: entry,
                             route: route, trans_logprob: tlog,
                             obs_index: obs_index)
          end
          entries << best if best
        end
        entries
      end

      def build_segment(seg_start, layer)
        final = layer.first
        layer.each { |e| final = e if e.score > final.score + SCORE_EPS }
        entries = []
        entry = final
        while entry
          entries.unshift(entry)
          entry = entry.back
        end
        matches = entries.map do |en|
          Match.new(en.obs_index, en.candidate.edge.id, en.candidate.fraction,
                    en.candidate.point, en.candidate.distance,
                    en.candidate.emission_logprob)
        end
        transitions = []
        entries.each_cons(2) do |a, b|
          transitions << Transition.new(a.obs_index, b.obs_index, b.route,
                                        b.trans_logprob)
        end
        emission = entries.sum { |en| en.candidate.emission_logprob }
        transition = transitions.sum(&:logprob)
        Segment.new(seg_start, entries.last.obs_index, matches, transitions,
                    final.score, emission, transition)
      end
    end
  end
end
