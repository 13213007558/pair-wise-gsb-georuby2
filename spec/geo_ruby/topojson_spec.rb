require File.expand_path('../../spec_helper', __FILE__)

describe GeoRuby::TopoJSONExporter do
  let(:sf) { GeoRuby::SimpleFeatures }

  def export(geometry, quantization = 3, name: 'collection')
    described_class.new(geometry, object_name: name,
                                  quantization: quantization).export
  end

  def decode_arc(index, topology)
    encoded = index >= 0 ? topology['arcs'][index] : topology['arcs'][~index]
    points = [encoded.first.dup]
    encoded[1..-1].each do |delta|
      points << points.last.each_with_index.map { |value, i| value + delta[i] }
    end
    points.reverse! if index.negative?
    points
  end

  def merge_arc_references(references, topology, ring: false)
    points = references.flat_map { |reference| decode_arc(reference, topology) }
    positions = []
    points.each do |point|
      positions << point if positions.empty? || positions.last != point
    end
    if ring && positions.any? && positions.last != positions.first
      positions << positions.first
    end
    positions
  end

  def decode_polygon(geometry, topology)
    geometry['arcs'].map do |ring|
      merge_arc_references(ring, topology, ring: true)
    end
  end

  let(:left_polygon) do
    sf::Polygon.from_coordinates([
      [[0, 0], [0, 2], [1, 2], [1, 0], [0, 0]]
    ])
  end

  let(:right_polygon) do
    sf::Polygon.from_coordinates([
      [[1, 0], [2, 0], [2, 2], [1, 2], [1, 0]]
    ])
  end

  it 'quantizes before topology and shares a common edge by one arc' do
    features = GeoRuby::GeoJSONFeatureCollection.new([
      GeoRuby::GeoJSONFeature.new(left_polygon, { 'name' => 'left' }, 7),
      GeoRuby::GeoJSONFeature.new(right_polygon, { 'name' => 'right' }, 8)
    ])
    result = export(features)
    topology = result.topology
    collection = topology['objects']['collection']
    left = collection['geometries'][0]
    right = collection['geometries'][1]

    expect(result.diagnostics).to eq([])
    expect(left['id']).to eq(7)
    expect(left['properties']).to eq('name' => 'left')
    expect(right['id']).to eq(8)
    expect(right['properties']).to eq('name' => 'right')

    shared_left = left['arcs'][0].find { |arc| arc.negative? }
    shared_right = right['arcs'][0].find { |arc| arc.abs == shared_left.abs }
    expect(shared_right).to eq(shared_left.abs)

    left_coordinates = decode_polygon(left, topology).first
    right_coordinates = decode_polygon(right, topology).first
    expected_left = [[0, 0], [0, 2], [1, 2], [1, 0], [0, 0]]
    expected_right = [[1, 0], [2, 0], [2, 2], [1, 2], [1, 0]]
    expect(left_coordinates).to eq([expected_left[2], expected_left[3],
                                   expected_left[0], expected_left[1], expected_left[2]])
    expect(right_coordinates).to eq([
      expected_right[1], expected_right[2], expected_right[3],
      expected_right[0], expected_right[1]
    ])
    expect(shared_left.abs).to eq(shared_right.abs)
    shared_arc = decode_arc(shared_left.abs, topology)
    expect(shared_arc).to eq([[1, 2], [1, 0]])
  end

  it 'matches boundaries with a rotated start and unequal input segmentation' do
    rotated_left = sf::Polygon.from_coordinates([
      [[1, 2], [1, 0], [0, 0], [0, 2], [1, 2]]
    ])
    segmented_right = sf::Polygon.from_coordinates([
      [[1, 0], [1, 1], [2, 1], [2, 2], [1, 2], [1, 1], [1, 0]]
    ])
    topology = export(GeoRuby::GeoJSONFeatureCollection.new([
      GeoRuby::GeoJSONFeature.new(rotated_left, {}, 'a'),
      GeoRuby::GeoJSONFeature.new(segmented_right, {}, 'b')
    ])).topology
    geometries = topology['objects']['collection']['geometries']
    shared = geometries.map { |geometry| geometry['arcs'][0].map(&:abs) }.reduce(:&)

    expect(shared).not_to be_empty
    decoded = geometries.map { |geometry| decode_polygon(geometry, topology).first }
    expect(decoded.first).to include([1, 1])
    expect(decoded.last).to include([1, 1])
    expect(decoded.first).to eq([[1, 1], [1, 0], [0, 0], [0, 2], [1, 2], [1, 1]])
    expect(decoded.last).to eq([[1, 1], [2, 1], [2, 2], [1, 2], [1, 1], [1, 0], [1, 1]])
  end

  it 'preserves point and empty geometry objects and holes' do
    point = sf::Point.from_x_y(0.5, 1.5)
    empty_line = sf::LineString.from_coordinates([])
    empty_points = sf::MultiPoint.from_coordinates([])
    empty_polygon = sf::Polygon.from_coordinates([])
    polygon = sf::Polygon.from_coordinates([
      [[0, 0], [0, 4], [4, 4], [4, 0], [0, 0]],
      [[1, 1], [3, 1], [3, 3], [1, 3], [1, 1]]
    ])
    features = [
      point, empty_line, empty_points, empty_polygon, polygon
    ].each_with_index.map do |geometry, index|
      GeoRuby::GeoJSONFeature.new(geometry, { 'index' => index }, "id-#{index}")
    end
    result = export(GeoRuby::GeoJSONFeatureCollection.new(features), 5)
    geometries = result.topology['objects']['collection']['geometries']

    expect(result.diagnostics).to eq([])
    expect(geometries[0]).to include(
      'type' => 'Point', 'coordinates' => [1, 2]
    )
    expect(geometries[1]).to include('type' => 'LineString', 'arcs' => [])
    expect(geometries[2]).to include('type' => 'MultiPoint', 'coordinates' => [])
    expect(geometries[3]).to include('type' => 'Polygon', 'arcs' => [])
    expect(geometries[4]['arcs'].length).to eq(2)
    expect(decode_polygon(geometries[4], result.topology)).to eq([
      [[0, 0], [0, 4], [4, 4], [4, 0], [0, 0]],
      [[1, 1], [3, 1], [3, 3], [1, 3], [1, 1]]
    ])
  end

  it 'returns diagnostics for objects not expressible after quantization' do
    collapsed_line = sf::LineString.from_coordinates([[0, 0], [0.1, 0.1]])
    collapsed_ring = sf::Polygon.from_coordinates([
      [[0, 0], [0.1, 0], [0.1, 0.1], [0, 0]]
    ])
    valid_point = sf::Point.from_x_y(2, 2)
    features = [
      collapsed_line, collapsed_ring, valid_point
    ].each_with_index.map do |geometry, index|
      GeoRuby::GeoJSONFeature.new(geometry, { 'index' => index }, "id-#{index}")
    end
    result = export(GeoRuby::GeoJSONFeatureCollection.new(features), 3)
    geometries = result.topology['objects']['collection']['geometries']

    expect(geometries[0]['type']).to be_nil
    expect(geometries[1]['type']).to be_nil
    expect(geometries[2]).to include('type' => 'Point', 'coordinates' => [2, 2])
    expect(geometries[0]).to include('id' => 'id-0', 'properties' => { 'index' => 0 })
    reasons = result.diagnostics.map { |diagnostic| diagnostic[:reason] }
    expect(reasons).to eq([:degenerate_line, :degenerate_ring])
    expect(result.diagnostics.map { |diagnostic| diagnostic[:id] }).to eq(%w[id-0 id-1])
  end

  it 'preserves Z coordinates through the quantized shared arc representation' do
    ring = lambda do |*coordinates|
      linear_ring = sf::LinearRing.new(4326, true)
      coordinates.each do |x, y, z|
        linear_ring << sf::Point.from_x_y_z(x, y, z)
      end
      linear_ring
    end
    left = sf::Polygon.new(4326, true)
    left << ring.call([0, 0, 2], [0, 4, 2], [2, 4, 2], [2, 0, 2], [0, 0, 2])
    right = sf::Polygon.new(4326, true)
    right << ring.call([2, 0, 2], [4, 0, 2], [4, 4, 2], [2, 4, 2], [2, 0, 2])
    topology = export(GeoRuby::GeoJSONFeatureCollection.new([
      GeoRuby::GeoJSONFeature.new(left, {}, 'a'),
      GeoRuby::GeoJSONFeature.new(right, {}, 'b')
    ]), 5).topology
    geometries = topology['objects']['collection']['geometries']
    left_shared = geometries.first['arcs'][0].find(&:negative?)
    shared_arc = topology['arcs'][left_shared.abs]
    expect(shared_arc[1..-1].all? { |delta| delta.length == 3 }).to eq(true)
    expect(decode_polygon(geometries.first, topology).first.map(&:last).uniq).to eq([2])
  end
end
