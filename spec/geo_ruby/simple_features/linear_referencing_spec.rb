require File.expand_path(File.dirname(__FILE__) + '/../../spec_helper')

describe GeoRuby::SimpleFeatures::LinearReferencing do
  Point = GeoRuby::SimpleFeatures::Point
  LineString = GeoRuby::SimpleFeatures::LineString
  MultiLineString = GeoRuby::SimpleFeatures::MultiLineString

  def point(x, y, z, m, srid = 123)
    Point.from_x_y_z_m(x, y, z, m, srid)
  end

  def sample_line
    LineString.from_points([
      point(0, 0, 10, 0),
      point(2, 0, 20, 2),
      point(4, 0, 30, 0),
      point(5, 0, 40, 0)
    ], 123, true, true)
  end

  it 'locates every occurrence of an interior measure by M interpolation' do
    result = sample_line.at_measure(1)

    expect(result.points.length).to eq(2)
    expect(result.points.map(&:point)).to contain_exactly(
      point(1, 0, 15, 1),
      point(3, 0, 25, 1)
    )
    expect(result.segments).to be_empty
  end

  it 'separates an isolated endpoint from a constant-M segment' do
    result = sample_line.at_measure(0)

    expect(result.points.map { |hit| [hit.point.x, hit.point.m] }).to eq([[0, 0]])
    expect(result.segments.length).to eq(1)
    expect(result.segments.first.line.to_coordinates).to eq([[4, 0, 30, 0], [5, 0, 40, 0]])
    expect(result.segments.first.m_direction).to eq(:constant)
  end

  it 'clips every traversal of an interval and retains interpolated Z and M' do
    result = sample_line.between_measures(0, 1)

    expect(result.slices.length).to eq(3)
    expect(result.slices.map(&:kind)).to eq(%i[line line line])
    expect(result.slices.map { |slice| slice.geometry.to_coordinates }).to eq([
      [[0, 0, 10, 0], [1, 0, 15, 1]],
      [[3, 0, 25, 1], [4, 0, 30, 0]],
      [[4, 0, 30, 0], [5, 0, 40, 0]]
    ])
    expect(result.slices.map(&:m_range)).to eq([[0, 1], [0, 1], [0, 0]])
    expect(result.slices.map(&:m_direction)).to eq(%i[ascending descending constant])
    expect(result.geometry).to be_instance_of(MultiLineString)
  end

  it 'honors interval endpoint openness' do
    result = sample_line.between_measures(0, 0, include_start: false, include_end: false)
    expect(result.slices).to be_empty
    expect(result.geometry).to be_nil

    result = sample_line.between_measures(0, 1, include_end: false)
    expect(result.slices.map(&:kind)).to eq(%i[line line line])
    expect(result.slices.map { |slice| [slice.start_inclusive, slice.end_inclusive] }).to eq([
      [true, false], [false, true], [true, true]
    ])
    expect(result.slices.select { |slice| slice.kind == :point }).to be_empty
  end

  it 'retains SRID, Z and direction for clipped sublines' do
    result = sample_line.between_measures(0, 2)
    line_slice = result.slices.find { |slice| slice.kind == :line && slice.start_segment_index.zero? }

    expect(line_slice.geometry.srid).to eq(123)
    expect(line_slice.geometry).to be_with_z
    expect(line_slice.geometry).to be_with_m
    expect(line_slice.m_direction).to eq(:ascending)
  end

  it 'reports missing M values without interpolating through them' do
    missing_point = Point.new(123, false, true)
    missing_point.set_x_y(0, 0)
    missing_point.m = nil
    line = LineString.from_points([
      missing_point,
      Point.from_x_y_m(2, 0, 2, 123),
      Point.from_x_y_m(4, 0, 4, 123)
    ], 123, false, true)
    result = line.between_measures(1, 3)

    expect(result.slices.map(&:geometry).map(&:to_coordinates)).to eq([[[2.0, 0.0, 2.0], [3.0, 0.0, 3.0]]])
    expect(result.issues).to include(
      have_attributes(line_index: 0, segment_index: nil, vertex_index: 0)
    )
  end

  it 'reports points without an enabled M ordinate but keeps a valid adjacent vertex' do
    line = LineString.from_points([
      Point.from_x_y(0, 0, 123),
      point(2, 0, 20, 2)
    ], 123, true, true)
    result = line.at_measure(2)

    expect(result.points.map(&:point)).to eq([point(2, 0, 20, 2)])
    expect(result.issues).to contain_exactly(
      have_attributes(line_index: 0, segment_index: nil, vertex_index: 0)
    )
  end

  it 'reports repeated vertices and preserves zero-width hits' do
    line = LineString.from_points([
      Point.from_x_y_m(0, 0, 0, 123),
      Point.from_x_y_m(0, 0, 0, 123),
      Point.from_x_y_m(2, 0, 1, 123)
    ], 123, false, true)
    locate = line.at_measure(0)
    range = line.between_measures(0, 0)

    expect(locate.issues).to contain_exactly(
      have_attributes(line_index: 0, segment_index: 0, vertex_index: 0)
    )
    expect(locate.points.first.duplicate_vertex).to eq(true)
    expect(range.slices.first.duplicate_vertex).to eq(true)
  end

  it 'keeps multi-line junction metadata for locates and clips' do
    first = LineString.from_points([
      Point.from_x_y_m(0, 0, 0, 123),
      Point.from_x_y_m(1, 0, 1, 123)
    ], 123, false, true)
    second = LineString.from_points([
      Point.from_x_y_m(1, 0, 1, 123),
      Point.from_x_y_m(2, 0, 2, 123)
    ], 123, false, true)
    multi = MultiLineString.from_line_strings([first, second])

    locate = multi.at_measure(1)
    range = multi.between_measures(0.5, 1.5)

    expect(multi).to be_with_m
    expect(locate.points.first.junction).to eq(true)
    expect(locate.junctions.map(&:line_indices)).to eq([[0, 1]])
    expect(range.slices.length).to eq(2)
    expect(range.slices.map(&:junction)).to eq([true, true])
  end
end
