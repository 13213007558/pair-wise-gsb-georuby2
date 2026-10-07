require File.expand_path(File.dirname(__FILE__) + '/../../spec_helper')

module JointSimplifierSpecHelper
  def polygon(*rings)
    GeoRuby::SimpleFeatures::Polygon.from_coordinates(rings)
  end

  def xy_list(ring)
    ring.points.map { |point| [point.x, point.y] }
  end

  # Extracts the open chain of +ring+ (array of [x, y]) going from
  # +from+ to +to+ along the ring direction.
  def chain_between(ring_xy, from, to)
    points = ring_xy[0..-2]
    start = points.index(from)
    stop = points.index(to)
    return nil unless start && stop
    chain = []
    i = start
    loop do
      chain << points[i]
      break if i == stop
      i = (i + 1) % points.size
    end
    chain
  end

  def point_segment_distance(p, a, b)
    dx = b[0] - a[0]
    dy = b[1] - a[1]
    length_sq = dx * dx + dy * dy
    return Math.hypot(p[0] - a[0], p[1] - a[1]) if length_sq.zero?
    t = ((p[0] - a[0]) * dx + (p[1] - a[1]) * dy) / length_sq
    t = [[t, 0.0].max, 1.0].min
    Math.hypot(p[0] - a[0] - t * dx, p[1] - a[1] - t * dy)
  end

  def segment_distance(a, b, c, d)
    return 0.0 if segments_intersect?(a, b, c, d)
    [point_segment_distance(a, c, d), point_segment_distance(b, c, d),
     point_segment_distance(c, a, b), point_segment_distance(d, a, b)].min
  end

  def segments_intersect?(a, b, c, d)
    o1 = cross2d(a, b, c)
    o2 = cross2d(a, b, d)
    o3 = cross2d(c, d, a)
    o4 = cross2d(c, d, b)
    return true if (o1 * o2).negative? && (o3 * o4).negative?
    (o1.zero? && on_seg?(a, b, c)) || (o2.zero? && on_seg?(a, b, d)) ||
      (o3.zero? && on_seg?(c, d, a)) || (o4.zero? && on_seg?(c, d, b))
  end

  def cross2d(a, b, p)
    (b[0] - a[0]) * (p[1] - a[1]) - (b[1] - a[1]) * (p[0] - a[0])
  end

  def on_seg?(a, b, p)
    cross2d(a, b, p).abs <= 1e-9 &&
      p[0] >= [a[0], b[0]].min - 1e-9 && p[0] <= [a[0], b[0]].max + 1e-9 &&
      p[1] >= [a[1], b[1]].min - 1e-9 && p[1] <= [a[1], b[1]].max + 1e-9
  end

  # Minimum distance between the boundaries of two rings (arrays of
  # [x, y] including the closing point).
  def boundary_distance(ring_a, ring_b)
    ring_a.each_cons(2).map do |a, b|
      ring_b.each_cons(2).map { |c, d| segment_distance(a, b, c, d) }.min
    end.min
  end

  def point_to_polyline_distance(p, polyline)
    polyline.each_cons(2).map { |a, b| point_segment_distance(p, a, b) }.min
  end
end

describe GeoRuby::SimpleFeatures::JointSimplifier do
  include JointSimplifierSpecHelper

  let(:shared_edge) { [[5, 0], [5.3, 2], [5.9, 4], [5.8, 6], [5.3, 8], [5, 10]] }

  let(:west) do
    polygon([[0, 0], [5, 0], [5.3, 2], [5.9, 4], [5.8, 6], [5.3, 8],
             [5, 10], [0, 10], [0, 0]])
  end

  let(:east) do
    polygon([[5, 0], [10, 0], [10, 10], [5, 10], [5.3, 8], [5.8, 6],
             [5.9, 4], [5.3, 2], [5, 0]])
  end

  describe 'shared boundaries' do
    it 'keeps shared arcs identical in every polygon' do
      simplified = described_class.simplify([west, east], 0.5)
      west_chain = chain_between(xy_list(simplified[0][0]), [5, 0], [5, 10])
      east_chain = chain_between(xy_list(simplified[1][0]), [5, 10], [5, 0])
      expect(west_chain).to eq([[5, 0], [5.9, 4], [5, 10]])
      expect(east_chain).to eq(west_chain.reverse)
    end

    it 'supports a shared edge resegmented into more collinear edges' do
      # east shares the same geometric boundary but with extra vertices
      # inserted in the middle of the west shared edge segments.
      east_resegmented = polygon(
        [[5, 0], [10, 0], [10, 10], [5, 10], [5.15, 9], [5.3, 8],
         [5.8, 6], [5.85, 5], [5.9, 4], [5.3, 2], [5.15, 1], [5, 0]]
      )
      simplifier = described_class.new([west, east_resegmented], 0.5)
      simplified = simplifier.simplify
      west_chain = chain_between(xy_list(simplified[0][0]), [5, 0], [5, 10])
      east_chain = chain_between(xy_list(simplified[1][0]), [5, 10], [5, 0])
      expect(west_chain).to eq(east_chain.reverse)
      # every original vertex of both segmentations stays within epsilon
      (shared_edge + [[5.15, 9], [5.85, 5], [5.15, 1]]).each do |vertex|
        expect(point_to_polyline_distance(vertex, west_chain)).to be <= 0.5
      end
    end

    it 'supports a straight shared edge split into different edge counts' do
      straight_west = polygon([[0, 0], [5, 0], [5, 10], [0, 10], [0, 0]])
      split_east = polygon([[5, 0], [10, 0], [10, 10], [5, 10],
                            [5, 7.5], [5, 2.5], [5, 0]])
      simplified = described_class.simplify([straight_west, split_east], 0.5)
      west_chain = chain_between(xy_list(simplified[0][0]), [5, 0], [5, 10])
      east_chain = chain_between(xy_list(simplified[1][0]), [5, 10], [5, 0])
      expect(west_chain).to eq([[5, 0], [5, 10]])
      expect(east_chain).to eq([[5, 10], [5, 0]])
    end

    it 'does not introduce gaps along the shared boundary' do
      simplified = described_class.simplify([west, east], 0.5)
      distance = boundary_distance(xy_list(simplified[0][0]), xy_list(simplified[1][0]))
      expect(distance).to eq(0.0)
    end
  end

  describe 'junctions' do
    # three regions meeting at the single node (5, 5)
    let(:region_a) do
      polygon([[0, 0], [5, 0], [5.3, 2.5], [5, 5], [2.5, 5.3], [0, 5], [0, 0]])
    end
    let(:region_b) do
      polygon([[5, 0], [10, 0], [10, 5], [7.5, 5.3], [5, 5], [5.3, 2.5], [5, 0]])
    end
    let(:region_c) do
      polygon([[0, 5], [2.5, 5.3], [5, 5], [7.5, 5.3], [10, 5],
               [10, 10], [0, 10], [0, 5]])
    end

    it 'never removes a node shared by three regions' do
      simplified = described_class.simplify([region_a, region_b, region_c], 0.35)
      simplified.each do |region|
        expect(xy_list(region[0])).to include([5, 5])
      end
    end

    it 'keeps all three shared edges consistent' do
      simplified = described_class.simplify([region_a, region_b, region_c], 0.35)
      rings = simplified.map { |region| xy_list(region[0]) }
      expect(chain_between(rings[0], [5, 0], [5, 5]))
        .to eq(chain_between(rings[1], [5, 5], [5, 0]).reverse)
      expect(chain_between(rings[0], [5, 5], [0, 5]))
        .to eq(chain_between(rings[2], [0, 5], [5, 5]).reverse)
      expect(chain_between(rings[1], [10, 5], [5, 5]))
        .to eq(chain_between(rings[2], [5, 5], [10, 5]).reverse)
    end
  end

  describe 'anchors' do
    it 'never removes user specified anchor vertices' do
      simplified = described_class.simplify([west, east], 0.5, anchors: [[5.3, 2]])
      expect(xy_list(simplified[0][0])).to include([5.3, 2])
      expect(xy_list(simplified[1][0])).to include([5.3, 2])
    end

    it 'accepts Point anchors and inserts anchors lying on a shared edge' do
      straight_west = polygon([[0, 0], [5, 0], [5, 10], [0, 10], [0, 0]])
      straight_east = polygon([[5, 0], [10, 0], [10, 10], [5, 10], [5, 0]])
      anchor = GeoRuby::SimpleFeatures::Point.from_x_y(5, 4)
      simplified = described_class.simplify([straight_west, straight_east], 1,
                                            anchors: [anchor])
      expect(xy_list(simplified[0][0])).to include([5, 4])
      expect(xy_list(simplified[1][0])).to include([5, 4])
    end

    it 'warns about anchors that are not on any boundary' do
      simplifier = described_class.new([west, east], 0.5, anchors: [[100, 100]])
      simplifier.simplify
      expect(simplifier.warnings.join).to include('anchor (100, 100) ignored')
    end
  end

  describe 'topology constraints' do
    # two faces separated by a narrow channel (0.3 wide at y = 5)
    let(:channel_west) do
      polygon([[0, 0], [10, 0], [8, 5], [10, 10], [0, 10], [0, 0]])
    end
    let(:channel_east) do
      polygon([[10.3, 0], [20, 0], [20, 10], [10.3, 10], [9.5, 5], [10.3, 0]])
    end

    it 'keeps chains whose simplification would touch the other face' do
      simplifier = described_class.new([channel_west, channel_east], 2.5)
      simplified = simplifier.simplify
      expect(simplifier.warnings).not_to be_empty
      expect(xy_list(simplified[0][0])).to include([8, 5])
    end

    it 'never makes originally separate faces touch' do
      simplified = described_class.simplify([channel_west, channel_east], 2.5)
      distance = boundary_distance(xy_list(simplified[0][0]), xy_list(simplified[1][0]))
      expect(distance).to be > 0
    end

    # a hole whose top vertex sits 0.15 below the outer boundary bulge
    let(:holed) do
      polygon(
        [[0, 0], [10, 0], [10, 10], [5, 10.2], [0, 10], [0, 0]],
        [[4, 9.9], [6, 9.9], [5, 10.05], [4, 9.9]]
      )
    end

    it 'keeps the outer boundary when simplification would swallow a nearby hole' do
      simplifier = described_class.new([holed], 0.5)
      simplified = simplifier.simplify
      expect(simplifier.warnings).not_to be_empty
      expect(xy_list(simplified[0][0])).to include([5, 10.2])
    end

    it 'keeps the hole inside its polygon and unchanged in membership' do
      simplified = described_class.simplify([holed], 0.5)
      outer = xy_list(simplified[0][0])
      hole = xy_list(simplified[0][1])
      expect(hole.uniq.size).to be >= 3
      expect(boundary_distance(outer, hole)).to be > 0
      hole[0..-2].each do |vertex|
        point = GeoRuby::SimpleFeatures::Point.from_x_y(vertex[0], vertex[1])
        expect(simplified[0][0].contains_point?(point)).to be(true)
      end
    end
  end

  describe 'result shape' do
    it 'returns polygons of the same class with preserved ring count' do
      simplified = described_class.simplify([west, east], 0.35)
      expect(simplified.size).to eq(2)
      simplified.each do |result|
        expect(result).to be_instance_of(GeoRuby::SimpleFeatures::Polygon)
        expect(result.rings.size).to eq(1)
        expect(result[0].closed?).to be(true)
      end
    end

    it 'accepts a MultiPolygon and returns one' do
      multi = GeoRuby::SimpleFeatures::MultiPolygon.from_polygons([west, east])
      simplified = described_class.simplify(multi, 0.35)
      expect(simplified).to be_instance_of(GeoRuby::SimpleFeatures::MultiPolygon)
      expect(simplified.geometries.size).to eq(2)
    end

    it 'reduces the total vertex count' do
      original = west[0].points.size + east[0].points.size
      simplified = described_class.simplify([west, east], 0.35)
      reduced = simplified.sum { |poly| poly[0].points.size }
      expect(reduced).to be < original
    end
  end
end
