require 'spec_helper'
require 'geo_ruby/map_matcher'

describe GeoRuby::MapMatcher do
  MM = GeoRuby::MapMatcher

  def network(edges)
    MM::RoadNetwork.new(edges)
  end

  def edge(id, points, oneway = false)
    { id: id, points: points, oneway: oneway }
  end

  describe 'candidate generation' do
    it 'projects onto segment interiors, not only vertices' do
      net = network([edge(:l, [[0, 0], [50, 0], [50, 50]])])
      frac, dist, pt = net.edges.first.project(50, 25)
      expect(frac).to be_within(1e-9).of(0.75)
      expect(dist).to be_within(1e-9).of(0.0)
      expect(pt).to eq([50.0, 25.0])

      cand = net.candidates(60, 25, 20, 5).first
      expect(cand.fraction).to be_within(1e-9).of(0.75)
      expect(cand.distance).to be_within(1e-9).of(10.0)
      expect(cand.point).to eq([50.0, 25.0])
    end

    it 'orders candidates deterministically and caps the count' do
      net = network([edge(:b, [[0, 0], [100, 0]]),
                     edge(:a, [[0, 3], [100, 3]]),
                     edge(:c, [[0, 6], [100, 6]])])
      cands = net.candidates(50, 0, 50, 2)
      expect(cands.map { |c| c.edge.id }).to eq(%i[b a])
    end
  end

  describe 'parallel roads' do
    # Two parallel two-way roads 20m apart, joined only by connectors at
    # x = 0 and x = 200, so switching roads mid-way is a long detour.
    let(:net) do
      network([edge(:a1, [[0, 0], [100, 0]]), edge(:a2, [[100, 0], [200, 0]]),
               edge(:b1, [[0, 20], [100, 20]]),
               edge(:b2, [[100, 20], [200, 20]]),
               edge(:c1, [[0, 0], [0, 20]]), edge(:c2, [[200, 0], [200, 20]])])
    end
    let(:obs) do
      [[10, 0], [30, 0], [50, 0], [70, 0], [90, 0],
       [110, 12], # noisy: closer to road B (8m) than to road A (12m)
       [130, 0], [150, 0], [170, 0], [190, 0]]
    end

    it 'stays on the driven road instead of jumping to the nearest one' do
      nearest = net.candidates(110, 12, 30, 8).first
      expect(nearest.edge.id).to eq(:b2) # point-wise nearest would jump

      result = MM::Matcher.new(net, sigma: 5, beta: 30).match(obs)
      expect(result.segments.length).to eq(1)
      seg = result.segments.first
      expect(seg.matches.map(&:obs_index)).to eq((0...obs.length).to_a)
      seg.matches.each do |m|
        expect(%i[a1 a2]).to include(m.edge_id)
      end
      seg.transitions.each do |t|
        expect(t.route.edge_ids - %i[a1 a2]).to be_empty
      end
    end
  end

  describe 'overpass without connection' do
    # road_b crosses road_a at (50, 0) but shares no vertex there.
    let(:net) do
      network([edge(:road_a, [[0, 0], [100, 0]]),
               edge(:road_b, [[50, -50], [50, 50]])])
    end

    it 'does not teleport onto the crossing road' do
      cands = net.candidates(50, 0, 10, 8)
      expect(cands.map { |c| c.edge.id }).to contain_exactly(:road_a, :road_b)

      obs = [[10, 0], [30, 0], [50, 0], [70, 0], [90, 0]]
      result = MM::Matcher.new(net, sigma: 5, beta: 30).match(obs)
      expect(result.segments.length).to eq(1)
      seg = result.segments.first
      expect(seg.matches.map(&:edge_id).uniq).to eq([:road_a])
      expect(seg.transitions.map { |t| t.route.length })
        .to all(be_within(1e-9).of(20.0))
    end
  end

  describe 'missing observation' do
    let(:net) { network([edge(:r, [[0, 0], [100, 0]], true)]) }

    it 'bridges the gap and reports the observation as unmatched' do
      obs = [{ t: 0, x: 10, y: 0 }, { t: 1, x: 30, y: 0 },
             { t: 2, x: nil, y: nil },
             { t: 3, x: 70, y: 0 }, { t: 4, x: 90, y: 0 }]
      result = MM::Matcher.new(net, sigma: 5, beta: 30).match(obs)
      expect(result.unmatched).to eq([MM::Unmatched.new(2, :missing)])
      expect(result.segments.length).to eq(1)
      seg = result.segments.first
      expect(seg.matches.map(&:obs_index)).to eq([0, 1, 3, 4])
      bridge = seg.transitions.find { |t| t.from_index == 1 }
      expect(bridge.to_index).to eq(3)
      expect(bridge.route.length).to be_within(1e-9).of(40.0)
    end

    it 'exposes per-segment scores that add up' do
      obs = [[10, 0], [30, 0], [70, 0]]
      result = MM::Matcher.new(net, sigma: 5, beta: 30).match(obs)
      seg = result.segments.first
      expect(seg.score).to be_within(1e-9)
        .of(seg.emission_score + seg.transition_score)
      expect(seg.transitions.map(&:logprob).sum)
        .to be_within(1e-9).of(seg.transition_score)
    end
  end

  describe 'disconnected components' do
    let(:net) do
      network([edge(:c1, [[0, 0], [100, 0]]),
               edge(:c2, [[0, 500], [100, 500]])])
    end

    it 'splits into explicit segments instead of drawing straight lines' do
      obs = [[10, 0], [30, 0], [30, 500], [50, 500], [30, 0]]
      result = MM::Matcher.new(net, sigma: 5, beta: 30).match(obs)
      expect(result.breaks).to eq([2, 4])
      expect(result.segments.length).to eq(3)
      expect(result.segments[0].matches.map(&:edge_id).uniq).to eq([:c1])
      expect(result.segments[1].matches.map(&:edge_id).uniq).to eq([:c2])
      expect(result.segments[2].matches.map(&:edge_id).uniq).to eq([:c1])
      # routes only ever exist inside a segment
      expect(result.segments[2].transitions).to be_empty
      result.segments.each { |s| expect(s.score).to be_finite }
    end

    it 'reports observations with no candidate in range' do
      obs = [[10, 0], [10, 250], [30, 0]]
      result = MM::Matcher.new(net, sigma: 5, beta: 30).match(obs)
      expect(result.unmatched).to eq([MM::Unmatched.new(1, :no_candidate)])
      expect(result.segments.length).to eq(2)
    end
  end

  describe 'direction rules' do
    it 'forbids reversing on a one-way edge without a network loop' do
      net = network([edge(:one, [[0, 0], [100, 0]], true)])
      result = MM::Matcher.new(net, sigma: 5, beta: 30)
                          .match([[80, 0], [20, 0]])
      expect(result.breaks).to eq([1])
      expect(result.segments.length).to eq(2)
    end

    it 'allows reversing on a two-way edge at travelled-distance cost' do
      net = network([edge(:two, [[0, 0], [100, 0]])])
      result = MM::Matcher.new(net, sigma: 5, beta: 30)
                          .match([[80, 0], [20, 0]])
      expect(result.segments.length).to eq(1)
      route = result.segments.first.transitions.first.route
      expect(route.length).to be_within(1e-9).of(60.0)
      step = route.steps.first
      expect(step.from_fraction).to be_within(1e-9).of(0.8)
      expect(step.to_fraction).to be_within(1e-9).of(0.2)
    end

    it 'adds a configurable u-turn penalty to same-edge reversal' do
      net = network([edge(:two, [[0, 0], [100, 0]])])
      matcher = MM::Matcher.new(net, sigma: 5, beta: 30, u_turn_penalty: 50)
      route = matcher.match([[80, 0], [20, 0]])
                     .segments.first.transitions.first.route
      expect(route.length).to be_within(1e-9).of(110.0)
    end

    it 'routes a one-way reversal around a network loop when one exists' do
      net = network([edge(:fwd, [[0, 0], [100, 0]], true),
                     edge(:back, [[100, 0], [100, 10], [0, 10], [0, 0]],
                          true)])
      result = MM::Matcher.new(net, sigma: 5, beta: 1000)
                          .match([[80, 0], [20, 0]])
      expect(result.segments.length).to eq(1)
      route = result.segments.first.transitions.first.route
      expect(route.steps.map(&:edge_id)).to eq(%i[fwd back fwd])
      expect(route.length).to be_within(1e-9).of(20 + 120 + 20)
    end
  end

  describe 'equivalent optima' do
    # Two identical roads on top of each other: every complete path
    # through either road scores exactly the same.
    let(:net) do
      network([edge(:b_edge, [[0, 0], [100, 0]]),
               edge(:a_edge, [[0, 0], [100, 0]])])
    end

    it 'breaks ties deterministically by candidate order' do
      matcher = MM::Matcher.new(net, sigma: 5, beta: 30)
      obs = [[10, 0], [20, 0], [30, 0]]
      first = matcher.match(obs)
      second = matcher.match(obs)
      expect(first.segments.length).to eq(1)
      ids = first.segments.first.matches.map(&:edge_id)
      expect(ids).to eq(%i[a_edge a_edge a_edge])
      expect(second.segments.first.matches.map(&:edge_id)).to eq(ids)
      expect(second.segments.first.score)
        .to eq(first.segments.first.score)
    end
  end
end
