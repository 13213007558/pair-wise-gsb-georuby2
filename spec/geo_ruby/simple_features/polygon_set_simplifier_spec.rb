require File.expand_path(File.dirname(__FILE__) + '/../../spec_helper')

describe GeoRuby::SimpleFeatures::PolygonSetSimplifier do
  def polygon(coordinates)
    GeoRuby::SimpleFeatures::Polygon.from_coordinates(coordinates)
  end

  def multi_polygon(*coordinates)
    GeoRuby::SimpleFeatures::MultiPolygon.from_polygons(coordinates.map { |rings| polygon(rings) })
  end

  it 'uses the same retained points on coincident reversed boundaries' do
    shared = [[0, 0], [1, 0.2], [2, 0], [3, 0.2], [4, 0]]
    bottom = shared + [[4, -4], [0, -4]]
    top = [[0, 0], [0, 4], [4, 4]] + shared.reverse
    geometry = multi_polygon([bottom], [top])

    simplified = geometry.simplify(0.3)
    bottom_boundary = simplified[0].rings.first.points.map(&:to_xy)
    top_boundary = simplified[1].rings.first.points.map(&:to_xy)

    expect(bottom_boundary).to include([0, 0], [2, 0], [4, 0])
    expect(top_boundary).to include([0, 0], [2, 0], [4, 0])
    expect(bottom_boundary).not_to include([1, 0.2], [3, 0.2])
    expect(top_boundary).not_to include([1, 0.2], [3, 0.2])

    shared_bottom = [[0, 0], [2, 0], [4, 0]]
    expect(bottom_boundary.each_cons(3).any? { |chain| chain == shared_bottom }).to be(true)
    expect(top_boundary.each_cons(3).any? { |chain| chain == shared_bottom.reverse }).to be(true)
  end

  it 'supports different collinear subdivisions on a shared edge' do
    left = polygon([[[0, 0], [0, 4], [4, 4], [4, 0], [0, 0]]])
    right = polygon([[[4, 0], [4, 2], [4, 4], [8, 4], [8, 0], [4, 0]]])
    geometry = GeoRuby::SimpleFeatures::MultiPolygon.from_polygons([left, right])

    simplified = geometry.simplify(2)

    expect(simplified[0].rings.first.points.map(&:to_xy).uniq)
      .to contain_exactly([0, 0], [0, 4], [4, 4], [4, 0])
    expect(simplified[1].rings.first.points.map(&:to_xy).uniq)
      .to contain_exactly([4, 0], [4, 4], [8, 4], [8, 0])
  end

  it 'keeps a node shared by three regions' do
    southwest = polygon([[[4, 4], [4, 0.2], [0, 0], [0, 8], [4, 4]]])
    southeast = polygon([[[4, 4], [8, 0], [8, 8], [4, 4]]])
    north = polygon([[[4, 4], [0, 8], [8, 8], [4, 4]]])
    geometry = GeoRuby::SimpleFeatures::MultiPolygon.from_polygons([southwest, southeast, north])

    simplified = geometry.simplify(1)
    all_coordinates = simplified.geometries.flat_map { |poly| poly.rings.flat_map { |ring| ring.points.map(&:to_xy) } }

    expect(all_coordinates).to include([4, 4], [4, 0.2])
  end

  it 'keeps user-specified anchors on a shared boundary' do
    shared = [[0, 0], [1, 0.2], [2, 0], [3, 0.2], [4, 0]]
    bottom = shared + [[4, -4], [0, -4]]
    top = [[0, 0], [0, 4], [4, 4]] + shared.reverse
    geometry = multi_polygon([bottom], [top])

    simplified = geometry.simplify(0.3, anchors: [GeoRuby::SimpleFeatures::Point.from_x_y(3, 0.2)])
    bottom_boundary = simplified[0].rings.first.points.map(&:to_xy)
    top_boundary = simplified[1].rings.first.points.map(&:to_xy)

    expect(bottom_boundary).to include([3, 0.2])
    expect(top_boundary).to include([3, 0.2])
  end

  it 'raises when an anchor is not on a polygon boundary' do
    geometry = multi_polygon([[[0, 0], [4, 0], [4, 4], [0, 4], [0, 0]]])

    expect { geometry.simplify(1, anchors: [[2, 5]]) }
      .to raise_error(ArgumentError, /not on a polygon boundary/)
  end

  it 'does not close a narrow channel' do
    outer = [[0, 0], [8, 0], [8, 10], [0, 10], [0, 6], [3.9, 6], [3.9, 4], [0, 4], [0, 0]]
    hole = [[5, 5], [7, 5], [7, 9], [5, 9], [5, 5]]
    geometry = multi_polygon([outer, hole])

    result = geometry.simplify(1.5, report: true)
    coordinates = result.multi_polygon[0].rings.first.points.map(&:to_xy)

    expect(coordinates).to include([3.9, 6], [3.9, 4], [0, 6], [0, 4])
    expect(result.warnings).not_to be_empty
    expect(result.warnings.first[:reason]).to eq('topology_constraint')
  end

  it 'does not let separate faces touch after simplification' do
    first = [[0, 0], [5, 0.9], [10, 0], [10, 4], [0, 4], [0, 0]]
    second = [[4, 5], [9, 5], [9, 9], [4, 9], [4, 5]]
    geometry = multi_polygon([first], [second])

    result = geometry.simplify(5, report: true)

    expect(result.multi_polygon[0].rings.first.points.map(&:to_xy)).to eq(first)
    expect(result.warnings.map { |warning| warning[:reason] }.uniq).to eq(['topology_constraint'])
  end

  it 'preserves a hole close to an outer boundary' do
    outer = [[0, 0], [0, 10], [10, 10], [10, 0], [5, 0.4], [0, 0]]
    hole = [[4, 2], [6, 2], [6, 4], [4, 4], [4, 2]]
    geometry = multi_polygon([outer, hole])

    result = geometry.simplify(1.5, report: true)

    expect(result.multi_polygon[0].rings.first.points.map(&:to_xy)).to include([5, 0.4])
    expect(result.multi_polygon[0].rings[1].points.map(&:to_xy)).to include([4, 2], [6, 2], [6, 4], [4, 4])
    expect(result.warnings).not_to be_empty
  end
end
