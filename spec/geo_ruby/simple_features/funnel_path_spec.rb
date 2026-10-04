require File.expand_path(File.dirname(__FILE__) + '/../../spec_helper')

describe GeoRuby::SimpleFeatures::FunnelPath do
  def point(x, y)
    GeoRuby::SimpleFeatures::Point.from_x_y(x.to_f, y.to_f)
  end

  def cross(first, second, third)
    (second.x - first.x) * (third.y - first.y) -
      (second.y - first.y) * (third.x - first.x)
  end

  def on_segment?(first, second, point)
    return false unless cross(first, second, point).zero?
    point.x >= [first.x, second.x].min && point.x <= [first.x, second.x].max &&
      point.y >= [first.y, second.y].min && point.y <= [first.y, second.y].max
  end

  def segments_intersect?(first, second, third, fourth)
    orientations = [
      cross(third, fourth, first),
      cross(third, fourth, second),
      cross(first, second, third),
      cross(first, second, fourth)
    ]

    return true if on_segment?(third, fourth, first)
    return true if on_segment?(third, fourth, second)
    return true if on_segment?(first, second, third)
    return true if on_segment?(first, second, fourth)

    orientations[0] * orientations[1] < 0 && orientations[2] * orientations[3] < 0
  end

  def event_indices(path, portals)
    path.points[1..-2].map do |path_point|
      portals.index do |left, right|
        path_point == left || path_point == right
      end
    end
  end

  def stay_in_corridor?(path, start_point, end_point, portals)
    indices = [-1] + event_indices(path, portals) + [portals.length]
    path.points.each_cons(2).with_index.all? do |(first, second), segment_index|
      (indices[segment_index] + 1...indices[segment_index + 1]).all? do |portal_index|
        left, right = portals[portal_index]
        next true if left == right
        segments_intersect?(first, second, left, right)
      end
    end
  end

  def visibility_graph_length(start_point, end_point, portals)
    nodes = [{ point: start_point, index: -1 }]
    portals.each_with_index do |(left, right), index|
      nodes << { point: left, index: index }
      nodes << { point: right, index: index } unless left == right
    end
    nodes << { point: end_point, index: portals.length }

    distances = Array.new(nodes.length) { Float::INFINITY }
    distances[0] = 0.0

    nodes.each_with_index do |source, source_index|
      nodes.each_with_index do |target, target_index|
        next unless target_index > source_index
        next unless target[:index] > source[:index]

        visible = (source[:index] + 1...target[:index]).all? do |portal_index|
          left, right = portals[portal_index]
          segments_intersect?(source[:point], target[:point], left, right)
        end
        next unless visible

        distance = distances[source_index] +
                   source[:point].euclidian_distance(target[:point])
        distances[target_index] = [distance, distances[target_index]].min
      end
    end

    distances.last
  end

  it 'builds a LineString through a straight portal corridor' do
    start_point = point(0, 0)
    end_point = point(10, 0)
    portals = [
      [point(2, -2), point(2, 2)],
      [point(2, -2), point(4, 2)],
      [point(4, -2), point(4, 2)],
      [point(4, 2), point(6, -2)],
      [point(6, -2), point(6, 2)],
      [point(6, 2), point(8, -2)]
    ]

    path = described_class.shortest_path(start_point, end_point, portals)

    expect(path).to be_a(GeoRuby::SimpleFeatures::LineString)
    expect(path.points).to eq([start_point, end_point])
    expect(path.euclidian_distance).to be_within(1e-10).of(10)
  end

  it 'matches the visibility graph shortest path through a folded corridor' do
    start_point = point(0, 0)
    end_point = point(10, 0)
    portals = [
      [point(2, -2), point(2, 2)],
      [point(4, 1), point(2, -2)],
      [point(6, 2), point(4, 1)],
      [point(8, -2), point(6, 2)]
    ]

    path = described_class.shortest_path(start_point, end_point, portals)

    expect(path.points.map(&:to_xy)).to eq([[0.0, 0.0], [4.0, 1.0], [10.0, 0.0]])
    expect(stay_in_corridor?(path, start_point, end_point, portals)).to be(true)
    expect(path.euclidian_distance)
      .to be_within(1e-10).of(visibility_graph_length(start_point, end_point, portals))
  end

  it 'normalizes arbitrary portal directions' do
    start_point = point(0, 0)
    end_point = point(10, 0)
    portals = [
      [point(2, 2), point(2, -2)],
      [point(2, -2), point(4, 1)],
      [point(4, 1), point(6, 2)],
      [point(8, -2), point(6, 2)]
    ]

    path = described_class.shortest_path(start_point, end_point, portals)

    expect(stay_in_corridor?(path, start_point, end_point, portals)).to be(true)
    expect(path.euclidian_distance)
      .to be_within(1e-10).of(visibility_graph_length(start_point, end_point, portals))
  end

  it 'forces the path through a zero-width portal' do
    start_point = point(0, 0)
    end_point = point(10, 0)
    pinch = point(2, -2)
    portals = [
      [point(2, 2), pinch],
      [pinch, pinch],
      [pinch, point(4, 2)]
    ]

    path = described_class.shortest_path(start_point, end_point, portals)

    expect(path.points).to eq([start_point, pinch, end_point])
    expect(stay_in_corridor?(path, start_point, end_point, portals)).to be(true)
    expect(path.euclidian_distance)
      .to be_within(1e-10).of(visibility_graph_length(start_point, end_point, portals))
  end

  it 'handles repeated and long collinear boundary portals' do
    start_point = point(0, 0)
    end_point = point(9, 0)
    first_portal = [point(1, 2), point(1, -2)]
    portals = Array.new(8) { [first_portal[0], first_portal[1]] }

    (1..7).each do |x|
      portals << [point(x, -2), point(x + 1, 2)]
      portals << [point(x + 1, 2), point(x + 1, -2)]
    end

    path = described_class.shortest_path(start_point, end_point, portals)

    expect(path.points).to eq([start_point, end_point])
    expect(stay_in_corridor?(path, start_point, end_point, portals)).to be(true)
    expect(path.euclidian_distance)
      .to be_within(1e-10).of(visibility_graph_length(start_point, end_point, portals))
  end

  it 'accepts coordinate arrays and exposes a LineString constructor' do
    path = GeoRuby::SimpleFeatures::LineString.funnel_path(
      [0, 0], [4, 0],
      [[[1, -1], [2, 1]], [[2, 1], [2, -1]]]
    )

    expect(path.points.map(&:to_xy)).to eq([[0.0, 0.0], [4.0, 0.0]])
  end

  it 'returns the direct segment for an empty portal list' do
    path = described_class.shortest_path(point(0, 0), point(3, 4), [])

    expect(path.points.map(&:to_xy)).to eq([[0.0, 0.0], [3.0, 4.0]])
    expect(path.euclidian_distance).to be_within(1e-10).of(5)
  end

  it 'forces zero-width portals at the beginning and end' do
    start_point = point(0, 0)
    end_point = point(4, 0)
    first_pin = point(1, 1)
    last_pin = point(3, -1)
    portals = [
      [first_pin, first_pin],
      [first_pin, point(2, -2)],
      [point(2, -2), point(3, 1)],
      [point(3, 1), last_pin],
      [last_pin, last_pin]
    ]

    path = described_class.shortest_path(start_point, end_point, portals)

    expect(path.points.map(&:to_xy))
      .to eq([[0.0, 0.0], [1.0, 1.0], [3.0, -1.0], [4.0, 0.0]])
    expect(stay_in_corridor?(path, start_point, end_point, portals)).to be(true)
  end

  it 'rejects portals without two endpoints' do
    expect do
      described_class.shortest_path(point(0, 0), point(1, 0), [[point(0, 0)]])
    end.to raise_error(described_class::Error, /two endpoints/)
  end

  it 'rejects non-contiguous portals' do
    portals = [
      [point(1, -1), point(1, 1)],
      [point(2, -1), point(2, 1)]
    ]

    expect do
      described_class.shortest_path(point(0, 0), point(3, 0), portals)
    end.to raise_error(described_class::Error, /do not share an endpoint/)
  end
end
