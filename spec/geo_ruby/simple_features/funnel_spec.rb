require File.expand_path('../../../spec_helper', __FILE__)

include GeoRuby::SimpleFeatures

describe Funnel do
  def pt(x, y)
    Point.from_x_y(x, y)
  end

  def coords(path)
    path.map { |p| [p.x, p.y] }
  end

  def path_length(path)
    path.each_cons(2).reduce(0.0) do |sum, (a, b)|
      sum + Math.hypot(b.x - a.x, b.y - a.y)
    end
  end

  # --- geometry helpers used to cross-check the funnel output ---

  def orient(a, b, c)
    (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x)
  end

  def on_segment?(a, b, p)
    orient(a, b, p) == 0 &&
      [a.x, b.x].min - 1e-12 <= p.x && p.x <= [a.x, b.x].max + 1e-12 &&
      [a.y, b.y].min - 1e-12 <= p.y && p.y <= [a.y, b.y].max + 1e-12
  end

  # true when the two segments cross at a point interior to both
  def proper_cross?(a, b, c, d)
    o1 = orient(a, b, c)
    o2 = orient(a, b, d)
    o3 = orient(c, d, a)
    o4 = orient(c, d, b)
    ((o1 > 0 && o2 < 0) || (o1 < 0 && o2 > 0)) &&
      ((o3 > 0 && o4 < 0) || (o3 < 0 && o4 > 0))
  end

  # closed polygon as an array of points (last == first)
  def corridor_polygon(start, goal, portals)
    ring = [start] + portals.map(&:first) + [goal] +
           portals.map(&:last).reverse
    deduped = ring.each_with_object([]) do |p, acc|
      acc << p if acc.empty? ||
                  acc.last.x != p.x || acc.last.y != p.y
    end
    deduped << deduped.first
    deduped
  end

  def point_in_polygon?(p, poly)
    return true if poly.each_cons(2).any? { |a, b| on_segment?(a, b, p) }

    inside = false
    poly.each_cons(2) do |a, b|
      next unless (a.y > p.y) != (b.y > p.y)

      x = a.x + (p.y - a.y) * (b.x - a.x) / (b.y - a.y)
      inside = !inside if x > p.x
    end
    inside
  end

  # a segment stays inside the (closed) corridor when no boundary edge
  # crosses it properly and sample points along it are inside
  def segment_in_corridor?(a, b, poly)
    [0.25, 0.5, 0.75].each do |t|
      sample = Point.from_x_y(a.x + t * (b.x - a.x),
                              a.y + t * (b.y - a.y))
      return false unless point_in_polygon?(sample, poly)
    end
    poly.each_cons(2) do |c, d|
      return false if proper_cross?(a, b, c, d)
    end
    true
  end

  # shortest path length through the corridor according to a
  # visibility graph over the corridor vertices (Dijkstra)
  def visibility_graph_length(start, goal, portals)
    poly = corridor_polygon(start, goal, portals)
    nodes = [start, goal] + portals.flatten(1)
    nodes = nodes.each_with_object([]) do |p, acc|
      acc << p unless acc.any? { |q| q.x == p.x && q.y == p.y }
    end
    n = nodes.length
    adj = Array.new(n) { [] }
    n.times do |i|
      n.times do |j|
        next if i == j
        next unless segment_in_corridor?(nodes[i], nodes[j], poly)

        w = Math.hypot(nodes[i].x - nodes[j].x, nodes[i].y - nodes[j].y)
        adj[i] << [j, w]
      end
    end
    dist = Array.new(n, Float::INFINITY)
    dist[0] = 0.0
    visited = Array.new(n, false)
    n.times do
      u = nil
      best = Float::INFINITY
      n.times do |k|
        next if visited[k] || dist[k] >= best

        best = dist[k]
        u = k
      end
      break if u.nil?

      visited[u] = true
      adj[u].each do |v, w|
        dist[v] = dist[u] + w if dist[u] + w < dist[v]
      end
    end
    dist[1]
  end

  # corridors used both for the visibility-graph comparison and for
  # the per-segment containment check
  def example_corridors
    [
      # straight strip, triangulated as a zig-zag: long runs of
      # collinear left and right endpoints
      {
        start: pt(-1, 0), goal: pt(9, 0),
        portals: [
          [pt(1, 1), pt(0, -1)],
          [pt(1, 1), pt(2, -1)],
          [pt(3, 1), pt(2, -1)],
          [pt(3, 1), pt(4, -1)],
          [pt(5, 1), pt(4, -1)],
          [pt(5, 1), pt(6, -1)],
          [pt(7, 1), pt(6, -1)],
          [pt(7, 1), pt(8, -1)]
        ]
      },
      # L-shaped corridor bending around the reflex vertex (2, 0)
      {
        start: pt(0, 1), goal: pt(3, -3),
        portals: [
          [pt(0, 2), pt(2, 0)],
          [pt(2, 2), pt(2, 0)],
          [pt(4, 2), pt(2, 0)],
          [pt(4, -3), pt(2, 0)]
        ]
      },
      # wavy corridor
      {
        start: pt(0, 0), goal: pt(7, 1),
        portals: [
          [pt(1, 2), pt(2, -1)],
          [pt(3, 2), pt(2, -1)],
          [pt(3, 2), pt(5, 0)],
          [pt(6, 3), pt(5, 0)]
        ]
      },
      # corridor pinched to a zero-width portal at (2, 0)
      {
        start: pt(0, 0), goal: pt(5, 0),
        portals: [
          [pt(1, 1), pt(1, -1)],
          [pt(1, 1), pt(2, 0)],
          [pt(2, 0), pt(2, 0)],
          [pt(2, 0), pt(3, -1)],
          [pt(4, 1), pt(3, -1)]
        ]
      },
      # S-shaped corridor with two bends
      {
        start: pt(0, 0), goal: pt(8, 1),
        portals: [
          [pt(1, 1), pt(1, -2)],
          [pt(3, 2), pt(1, -2)],
          [pt(3, 2), pt(4, -1)],
          [pt(5, 3), pt(4, -1)],
          [pt(5, 3), pt(6, 0)],
          [pt(7, 2), pt(6, 0)]
        ]
      }
    ]
  end

  describe 'portal validation' do
    it 'rejects discontinuous portals' do
      portals = [[pt(1, 1), pt(1, -1)], [pt(3, 1), pt(3, -1)]]
      expect { Funnel.path(pt(0, 0), pt(4, 0), portals) }
        .to raise_error(ArgumentError, /continuous/)
    end

    it 'rejects portals whose shared endpoint switches side' do
      portals = [[pt(1, 1), pt(1, -1)], [pt(1, -1), pt(3, 0)]]
      expect { Funnel.path(pt(0, 0), pt(4, 0), portals) }
        .to raise_error(ArgumentError, /orientation/)
    end

    it 'accepts zero-width portals' do
      portals = [
        [pt(1, 1), pt(2, 0)],
        [pt(2, 0), pt(2, 0)],
        [pt(2, 0), pt(3, -1)]
      ]
      expect { Funnel.new(portals) }.not_to raise_error
    end

    it 'accepts [x, y] arrays as points' do
      path = Funnel.path([0, 0], [4, 0],
                         [[[1, 1], [1, -1]], [[3, 1], [1, -1]]])
      expect(path.first).to be_a_point(0, 0)
      expect(path.last).to be_a_point(4, 0)
    end
  end

  describe '#path' do
    it 'returns the straight segment when there are no portals' do
      expect(coords(Funnel.path(pt(0, 0), pt(3, 2), [])))
        .to eq([[0, 0], [3, 2]])
    end

    it 'keeps a straight line through a long chain of collinear portals' do
      corridor = example_corridors[0]
      path = Funnel.path(corridor[:start], corridor[:goal],
                         corridor[:portals])
      expect(coords(path)).to eq([[-1, 0], [9, 0]])
    end

    it 'bends around the reflex vertex of an L-shaped corridor' do
      corridor = example_corridors[1]
      path = Funnel.path(corridor[:start], corridor[:goal],
                         corridor[:portals])
      expect(coords(path)).to eq([[0, 1], [2, 0], [3, -3]])
    end

    it 'passes through a zero-width portal' do
      corridor = example_corridors[3]
      path = Funnel.path(corridor[:start], corridor[:goal],
                         corridor[:portals])
      expect(coords(path)).to eq([[0, 0], [2, 0], [5, 0]])
    end

    it 'starts at start and ends at goal' do
      example_corridors.each do |corridor|
        path = Funnel.path(corridor[:start], corridor[:goal],
                           corridor[:portals])
        expect(path.first).to be_a_point(corridor[:start].x,
                                         corridor[:start].y)
        expect(path.last).to be_a_point(corridor[:goal].x,
                                        corridor[:goal].y)
      end
    end
  end

  describe 'cross-check against a visibility graph' do
    it 'matches the visibility-graph shortest path length' do
      example_corridors.each do |corridor|
        path = Funnel.path(corridor[:start], corridor[:goal],
                           corridor[:portals])
        expected = visibility_graph_length(corridor[:start],
                                           corridor[:goal],
                                           corridor[:portals])
        expect(path_length(path)).to be_within(1e-9).of(expected)
      end
    end

    it 'keeps every segment of the path inside the corridor' do
      example_corridors.each do |corridor|
        path = Funnel.path(corridor[:start], corridor[:goal],
                           corridor[:portals])
        poly = corridor_polygon(corridor[:start], corridor[:goal],
                                corridor[:portals])
        path.each_cons(2) do |a, b|
          expect(segment_in_corridor?(a, b, poly)).to be(true),
            "segment (#{a.x},#{a.y}) -> (#{b.x},#{b.y}) leaves the corridor"
        end
      end
    end
  end
end
