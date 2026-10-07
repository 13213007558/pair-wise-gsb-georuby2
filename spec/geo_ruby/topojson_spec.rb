require File.expand_path(File.dirname(__FILE__) + '/../spec_helper')
require 'geo_ruby/topojson'
require 'geo_ruby/topojson_decoder'

describe GeoRuby::TopoJSONExporter do
  let(:left_polygon) do
    GeoRuby::SimpleFeatures::Polygon.from_coordinates([[
      [0, 0], [0, 10], [10, 10], [10, 0], [0, 0]
    ]])
  end

  let(:right_polygon) do
    GeoRuby::SimpleFeatures::Polygon.from_coordinates([[
      [10, 10], [20, 10], [20, 0], [10, 0], [10, 10]
    ]])
  end


  let(:features) do
    [
      GeoRuby::GeoJSONFeature.new(left_polygon, { name: 'left' }, 'left-id'),
      GeoRuby::GeoJSONFeature.new(right_polygon, { name: 'right' }, 'right-id')
    ]
  end

  let(:shared_topology) { described_class.export(features, extent: 21) }

  it 'exports topology metadata, objects, properties and identity' do
    expect(shared_topology['type']).to eql('Topology')
    expect(shared_topology['transform']).to eql(
      'scale' => [1.0, 0.5],
      'translate' => [0.0, 0.0]
    )
    expect(shared_topology['arcs']).to be_an(Array)
    objects = shared_topology['objects']['data']['geometries']
    expect(objects.map { |object| object['id'] }).to eql(['left-id', 'right-id'])
    expect(objects.map { |object| object['properties']['name'] }).to eql(['left', 'right'])
  end

  it 'shares reverse arcs even when ring direction and start differ' do
    left_ring = shared_topology['objects']['data']['geometries'][0]['arcs'][0]
    right_ring = shared_topology['objects']['data']['geometries'][1]['arcs'][0]

    expect(left_ring).to include(-4)
    expect(right_ring).to include(2)
    expect(shared_topology['arcs'][2]).to eql([[10, 0], [0, 20]])
  end

  it 'shares arcs when one boundary is supplied in smaller segments' do
    segmented = GeoRuby::SimpleFeatures::LineString.from_coordinates([[0, 0], [0, 5], [0, 10]])
    continuous = GeoRuby::SimpleFeatures::LineString.from_coordinates([[0, 10], [0, 0]])
    segmented_topology = described_class.export([
      GeoRuby::GeoJSONFeature.new(segmented, {}),
      GeoRuby::GeoJSONFeature.new(continuous, {})
    ], extent: 11)
    objects = segmented_topology['objects']['data']['geometries']

    expect(objects[0]['arcs']).to eql([[0, 1]])
    expect(objects[1]['arcs']).to eql([[-2, -1]])
    expect(segmented_topology['arcs'].length).to eql(2)
  end

  it 'decodes delta arcs and closes rings to the agreed quantized coordinates' do
    decoded = GeoRuby::TopoJSONDecoder.coordinates(shared_topology, quantized: true)

    expect(decoded).to eql([
      [[[0, 0], [0, 20], [10, 20], [10, 0], [0, 0]]],
      [[[10, 20], [20, 20], [20, 0], [10, 0], [10, 20]]]
    ])
  end

  it 'preserves points, multipoints, holes and empty geometries' do
    point = GeoRuby::SimpleFeatures::Point.from_x_y(5, 2.5)
    multi_point = GeoRuby::SimpleFeatures::MultiPoint.from_coordinates([[0, 2.5], [10, 5.0]])
    geometry = GeoRuby::SimpleFeatures::GeometryCollection.from_geometries([
      point,
      multi_point,
      GeoRuby::SimpleFeatures::Polygon.from_coordinates([
        [[0, 0], [0, 10], [10, 10], [10, 0], [0, 0]],
        [[2, 2], [8, 2], [8, 8], [2, 8], [2, 2]]
      ]),
      GeoRuby::SimpleFeatures::LineString.from_coordinates([])
    ])
    geometry_topology = described_class.export(geometry, extent: 21)
    object = geometry_topology['objects']['data']

    expect(object['type']).to eql('GeometryCollection')
    point_transform = { scale: [1.0, 0.5], translate: [0.0, 0.0] }
    point_topology = described_class.export(point, transform: point_transform, extent: 21)['objects']['data']
    multi_point_transform = { scale: [1.0, 0.5], translate: [0.0, 0.0] }
    multi_point_topology = described_class.export(multi_point, transform: multi_point_transform, extent: 21)['objects']['data']
    expect(point_topology).to include('type' => 'Point', 'coordinates' => [5, 5])
    expect(multi_point_topology).to include('coordinates' => [[0, 5], [10, 10]])
    expect(object['geometries'][2]['arcs'].length).to eql(2)
    expect(object['geometries'][3]).to include('type' => 'LineString', 'arcs' => [])
  end

  it 'reports geometries that collapse under quantization' do
    line = GeoRuby::SimpleFeatures::LineString.from_coordinates([[0, 0], [0.4, 0.4]])
    transform = { scale: [1, 1], translate: [0, 0] }
    topology = described_class.export(line, transform: transform, extent: 11)

    expect(topology['objects']['data']).to include('type' => nil)
    expect(topology['diagnostics']).to include(
      hash_including('code' => 'degenerate_geometry')
    )
    expect(topology['arcs']).to eql([])
  end

  it 'reports coordinates outside an explicit transform extent' do
    transform = { scale: [1, 1], translate: [0, 0] }
    topology = described_class.export(GeoRuby::SimpleFeatures::Point.from_x_y(12, 0), transform: transform, extent: 11)

    expect(topology['objects']['data']).to include('type' => nil)
    expect(topology['diagnostics']).to include(hash_including('code' => 'coordinate_outside_extent'))
  end
end
