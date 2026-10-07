require File.expand_path(File.dirname(__FILE__) + '/../spec_helper')

module MapMatchingSpecHelper
  MM = GeoRuby::MapMatching

  def build_edge(id, coords, oneway: false, nodes: nil)
    MM::Edge.new(id, coords, oneway: oneway, nodes: nodes)
  end

  # Two parallel two-way roads A (y=0) and B (y=10), each split into ten
  # 100m edges, joined by a connector at each end. 22 edges total.
  def parallel_network
    edges = []
    10.times do |i|
      edges << build_edge("A#{i}", [[i * 100.0, 0.0], [(i + 1) * 100.0, 0.0]])
      edges << build_edge("B#{i}", [[i * 100.0, 10.0], [(i + 1) * 100.0, 10.0]])
    end
    edges << build_edge('C0', [[0.0, 0.0], [0.0, 10.0]])
    edges << build_edge('C1', [[1000.0, 0.0], [1000.0, 10.0]])
    MM::Network.new(edges)
  end

  # Ground road X along y=0 crossed at (500, 0) by overpass Y along x=500.
  # Y uses explicit node ids, so the crossing does not connect. 12 edges.
  def overpass_network
    edges = []
    10.times do |i|
      edges << build_edge("X#{i}", [[i * 100.0, 0.0], [(i + 1) * 100.0, 0.0]])
    end
    edges << build_edge('Y0', [[500.0, -500.0], [500.0, 0.0]], nodes: %i[y_s y_m])
    edges << build_edge('Y1', [[500.0, 0.0], [500.0, 500.0]], nodes: %i[y_m y_n])
    MM::Network.new(edges)
  end
end

MM = GeoRuby::MapMatching

describe GeoRuby::MapMatching do
  include MapMatchingSpecHelper

  describe GeoRuby::MapMatching::Edge do
    let(:edge) { build_edge('E', [[0.0, 0.0], [100.0, 0.0], [100.0, 100.0]]) }

    it 'projects on segment interiors, not only vertices' do
      proj = edge.project([50.0, 30.0])
      expect(proj[:offset]).to be_within(1e-9).of(50.0)
      expect(proj[:point]).to eq([50.0, 0.0])
      expect(proj[:distance]).to be_within(1e-9).of(30.0)
    end

    it 'projects across vertices of a multi-segment polyline' do
      proj = edge.project([110.0, 50.0])
      expect(proj[:offset]).to be_within(1e-9).of(150.0)
      expect(proj[:point]).to eq([100.0, 50.0])
      expect(proj[:distance]).to be_within(1e-9).of(10.0)
    end

    it 'clamps projections to the endpoints' do
      expect(edge.project([100.0, 150.0])[:offset]).to be_within(1e-9).of(200.0)
      expect(edge.project([-20.0, 0.0])[:offset]).to be_within(1e-9).of(0.0)
    end

    it 'interpolates points at offsets' do
      expect(edge.point_at(150.0)).to eq([100.0, 50.0])
      expect(edge.point_at(0.0)).to eq([0.0, 0.0])
    end
  end

  describe GeoRuby::MapMatching::Network do
    it 'finds candidates within the radius, nearest first' do
      network = parallel_network
      candidates = network.candidates([150.0, 8.0], 30.0)
      expect(candidates.map(&:edge_id)).to include('A1', 'B1')
      expect(candidates.first.edge_id).to eq('B1')
      expect(network.candidates([150.0, 100.0], 30.0)).to be_empty
    end

    context 'same-edge movement' do
      let(:network) { MM::Network.new([build_edge('T', [[0.0, 0.0], [1000.0, 0.0]])]) }

      def candidate_at(x)
        network.candidates([x, 0.0], 1.0).first
      end

      it 'routes forward along a two-way edge' do
        route = network.route(candidate_at(100.0), candidate_at(300.0))
        expect(route.length).to be_within(1e-9).of(200.0)
        expect(route.legs.size).to eq(1)
      end

      it 'allows reverse movement on a two-way edge' do
        route = network.route(candidate_at(300.0), candidate_at(100.0))
        expect(route.length).to be_within(1e-9).of(200.0)
        expect(route.legs.first.from_offset).to be > route.legs.first.to_offset
      end
    end

    context 'one-way edges' do
      let(:oneway) { build_edge('E', [[0.0, 0.0], [1000.0, 0.0]], oneway: true) }

      def candidate_at(network, x)
        network.candidates([x, 0.0], 1.0).first
      end

      it 'forbids reverse movement on a one-way edge without a loop' do
        network = MM::Network.new([oneway])
        expect(network.route(candidate_at(network, 300.0), candidate_at(network, 100.0))).to be_nil
      end

      it 'routes reverse movement around a network loop' do
        loop_edges = [
          oneway,
          build_edge('F', [[1000.0, 0.0], [1000.0, 100.0]], oneway: true),
          build_edge('G', [[1000.0, 100.0], [0.0, 100.0]], oneway: true),
          build_edge('H', [[0.0, 100.0], [0.0, 0.0]], oneway: true)
        ]
        network = MM::Network.new(loop_edges)
        route = network.route(candidate_at(network, 300.0), candidate_at(network, 100.0))
        expect(route.length).to be_within(1e-9).of(2000.0)
        expect(route.edge_ids).to eq(%w[E F G H E])
      end

      it 'allows forward movement on a one-way edge' do
        network = MM::Network.new([oneway])
        route = network.route(candidate_at(network, 100.0), candidate_at(network, 300.0))
        expect(route.length).to be_within(1e-9).of(200.0)
      end
    end

    context 'disconnected regions' do
      it 'returns nil when no route exists' do
        network = MM::Network.new([
                                    build_edge('P', [[0.0, 0.0], [500.0, 0.0]]),
                                    build_edge('Q', [[0.0, 10_000.0], [500.0, 10_000.0]])
                                  ])
        from = network.candidates([100.0, 0.0], 1.0).first
        to = network.candidates([100.0, 10_000.0], 1.0).first
        expect(network.route(from, to)).to be_nil
      end
    end

    context 'overpass without connection' do
      it 'sees the overpass as a candidate but cannot route to it' do
        network = overpass_network
        candidates = network.candidates([500.0, 2.0], 30.0)
        expect(candidates.map(&:edge_id)).to include('Y1')
        from = network.candidates([450.0, 2.0], 30.0).find { |c| c.edge_id == 'X4' }
        to = candidates.find { |c| c.edge_id == 'Y1' }
        expect(network.route(from, to)).to be_nil
      end
    end
  end

  describe GeoRuby::MapMatching::Matcher do
    context 'with noisy points bouncing between two parallel roads' do
      let(:network) { parallel_network }
      let(:matcher) { MM::Matcher.new(network, sigma: 5.0, beta: 30.0, max_radius: 30.0) }
      # True path runs on road A; noise pulls several fixes closer to road B.
      let(:noise) { [1.0, 8.0, 2.0, 9.0, 1.0, 7.0, 2.0, 8.0, 1.0, 9.0] }
      let(:observations) do
        noise.each_with_index.map { |y, i| { time: i * 10, point: [(i * 100.0) + 50.0, y] } }
      end
      let(:result) { matcher.match(observations) }

      it 'keeps the whole sequence on road A in one segment' do
        expect(result.unmatched).to be_empty
        expect(result.segments.size).to eq(1)
        segment = result.segments.first
        expect(segment.matched_points.map(&:edge_id)).to eq((0..9).map { |i| "A#{i}" })
        expect(segment.matched_points.map(&:offset)).to all(be_within(1e-9).of(50.0))
      end

      it 'would have jumped roads with per-point nearest matching' do
        expect(network.candidates([150.0, 8.0], 30.0).first.edge_id).to eq('B1')
        expect(result.segments.first.matched_points[1].edge_id).to eq('A1')
      end

      it 'returns routes between consecutive matched positions' do
        segment = result.segments.first
        expect(segment.routes.size).to eq(9)
        expect(segment.routes.map(&:length)).to all(be_within(1e-9).of(100.0))
        expect(segment.routes.first.edge_ids).to eq(%w[A0 A1])
      end

      it 'scores the segment with emissions and transitions' do
        segment = result.segments.first
        expected_emission = noise.sum { |y| (y**2) / (2.0 * 25.0) }
        expected_transition = noise.each_cons(2).sum do |a, b|
          (Math.hypot(100.0, b - a) - 100.0).abs / 30.0
        end
        expect(segment.transition_cost).to be_within(1e-9).of(expected_transition)
        expect(segment.emission_cost).to be_within(1e-9).of(expected_emission)
        expect(segment.score).to be_within(1e-9).of(expected_emission + expected_transition)
      end
    end

    context 'with an overpass crossing the driven road' do
      let(:network) { overpass_network }
      let(:matcher) { MM::Matcher.new(network, sigma: 5.0, beta: 30.0, max_radius: 30.0) }
      let(:observations) do
        (0..9).map { |i| { time: i * 10, point: [(i * 100.0) + 50.0, 2.0] } }
      end

      it 'never jumps onto the unconnected overpass' do
        result = matcher.match(observations)
        expect(result.segments.size).to eq(1)
        expect(result.segments.first.matched_points.map(&:edge_id)).to eq((0..9).map { |i| "X#{i}" })
      end
    end

    context 'with a missing observation' do
      let(:network) { MM::Network.new([build_edge('R', [[0.0, 0.0], [1000.0, 0.0]])]) }
      let(:matcher) { MM::Matcher.new(network, sigma: 5.0, beta: 30.0, max_radius: 30.0) }
      let(:observations) do
        (1..5).map { |i| { time: i * 10, point: i == 3 ? nil : [i * 100.0, 1.0] } }
      end
      let(:result) { matcher.match(observations) }

      it 'reports the gap and bridges it in a single segment' do
        expect(result.unmatched).to eq([2])
        expect(result.segments.size).to eq(1)
        segment = result.segments.first
        expect(segment.matched_points.size).to eq(4)
        expect(segment.matched_points.map(&:observation_index)).to eq([0, 1, 3, 4])
      end

      it 'routes across the gap using the surrounding observations' do
        segment = result.segments.first
        expect(segment.routes.map(&:length)).to eq([100.0, 200.0, 100.0])
        expect(segment.transition_cost).to be_within(1e-9).of(0.0)
      end
    end

    context 'with an observation too far from any road' do
      let(:network) { MM::Network.new([build_edge('R', [[0.0, 0.0], [1000.0, 0.0]])]) }
      let(:matcher) { MM::Matcher.new(network, sigma: 5.0, beta: 30.0, max_radius: 30.0) }

      it 'treats it as unmatched like a missing fix' do
        result = matcher.match([{ time: 0, point: [100.0, 1.0] },
                                { time: 10, point: [5000.0, 5000.0] },
                                { time: 20, point: [200.0, 1.0] }])
        expect(result.unmatched).to eq([1])
        expect(result.segments.size).to eq(1)
        expect(result.segments.first.matched_points.size).to eq(2)
      end
    end

    context 'with observations in two disconnected regions' do
      let(:network) do
        MM::Network.new([
                          build_edge('P', [[0.0, 0.0], [500.0, 0.0]]),
                          build_edge('Q', [[0.0, 10_000.0], [500.0, 10_000.0]])
                        ])
      end
      let(:matcher) { MM::Matcher.new(network, sigma: 5.0, beta: 30.0, max_radius: 30.0) }
      let(:observations) do
        [[100.0, 1.0], [200.0, 1.0], [300.0, 1.0],
         [100.0, 10_001.0], [200.0, 10_001.0], [300.0, 10_001.0]]
          .each_with_index.map { |point, i| { time: i * 10, point: point } }
      end
      let(:result) { matcher.match(observations) }

      it 'splits into explicit segments instead of inventing a straight link' do
        expect(result.segments.size).to eq(2)
        first, second = result.segments
        expect(first.matched_points.map(&:edge_id)).to eq(%w[P P P])
        expect(second.matched_points.map(&:edge_id)).to eq(%w[Q Q Q])
        expect(first.observation_range).to eq(0..2)
        expect(second.observation_range).to eq(3..5)
        expect(first.routes.map(&:edge_ids)).to eq([['P'], ['P']])
        expect(second.routes.map(&:edge_ids)).to eq([['Q'], ['Q']])
      end

      it 'scores each segment independently' do
        first, second = result.segments
        expected = 3 * (1.0 / (2.0 * 25.0))
        expect(first.score).to be_within(1e-9).of(expected)
        expect(second.score).to be_within(1e-9).of(expected)
        expect(first.transition_cost).to be_within(1e-9).of(0.0)
      end
    end

    context 'with equivalent optimal candidates' do
      let(:network) do
        MM::Network.new([
                          build_edge('P', [[0.0, 0.0], [100.0, 0.0]]),
                          build_edge('Q', [[0.0, 20.0], [100.0, 20.0]])
                        ])
      end
      let(:matcher) { MM::Matcher.new(network, sigma: 5.0, beta: 30.0, max_radius: 30.0) }

      it 'breaks ties deterministically by distance, edge id and offset' do
        first = matcher.match([{ time: 0, point: [50.0, 10.0] }])
        second = matcher.match([{ time: 0, point: [50.0, 10.0] }])
        expect(first.matched_points.first.edge_id).to eq('P')
        expect(second.matched_points.first.edge_id).to eq('P')
        expect(first.segments.first.score).to be_within(1e-9).of(100.0 / (2.0 * 25.0))
      end
    end
  end
end
