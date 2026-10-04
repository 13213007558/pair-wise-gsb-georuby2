require File.expand_path(File.dirname(__FILE__) + '/../../spec_helper')

describe GeoRuby::SimpleFeatures::RegularGridContourer do
  let(:contourer_class) { GeoRuby::SimpleFeatures::RegularGridContourer }

  def coordinates(result)
    result.to_a.map { |geometry| geometry.map(&:to_xy) }
  end

  it 'builds one line through intersections shared by adjacent cells' do
    values = [
      [-1, -1, 0.5],
      [1, 1, 2.5]
    ]
    result = contourer_class.new(values).contour(0.25)

    expect(result.line_strings.length).to eq(1)
    expect(result.linear_rings).to be_empty
    builder = contourer_class::LevelBuilder.allocate
    builder.send(:initialize, contourer_class.new(values), 0.25)
    builder.send(:add_cell, 0, 0)
    left_intersection = builder.point_cache[[:v, 0, 1, 0.625]]
    builder.send(:add_cell, 0, 1)
    right_intersection = builder.point_cache[[:v, 0, 1, 0.625]]
    expect(right_intersection).to equal(left_intersection)
    expect(coordinates(result)).to eq([
      [[0, 0.625], [1, 0.625], [1.8333333333333335, 0]]
    ])
  end

  it 'builds a closed ring around an interior high value' do
    values = [
      [-1, -1, -1],
      [-1, 1, -1],
      [-1, -1, -1]
    ]
    result = contourer_class.new(values).contour(0.5)

    expect(result.line_strings).to be_empty
    expect(result.linear_rings.length).to eq(1)
    ring = result.linear_rings.first
    expect(ring).to be_closed
    expect(ring.length).to eq(6)
    expect(ring.map(&:to_xy)).to eq([
      [1, 0.75], [0.75, 1], [1, 1.25], [1.25, 1], [1, 0.75], [1, 0.75]
    ])
  end

  it 'uses the asymptotic decider differently on each side of the saddle' do
    values = [
      [-3, 5],
      [4, -4]
    ]
    below = contourer_class.new(values).contour(0.4)
    at_saddle = contourer_class.new(values).contour(0.5)
    above = contourer_class.new(values).contour(0.6)

    expect(coordinates(below)).to eq([
      [[0.425, 0], [0, 0.4857142857142857]],
      [[1, 0.5111111111111111], [0.45, 1]]
    ])
    expect(coordinates(at_saddle)).to eq([
      [[0.4375, 0], [1, 0.5]],
      [[0, 0.5], [0.4375, 1]]
    ])
    expect(coordinates(above)).to eq([
      [[0.45, 0], [1, 0.48888888888888893]],
      [[0, 0.5142857142857143], [0.425, 1]]
    ])
  end

  it 'splits an exact saddle level into two crossing lines' do
    values = [
      [-1, 1],
      [1, -1]
    ]
    result = contourer_class.new(values).contour(0)

    expect(result.linear_rings).to be_empty
    expect(coordinates(result)).to contain_exactly(
      [[0.5, 0], [0.5, 0.5], [0.5, 1]],
      [[0, 0.5], [0.5, 0.5], [1, 0.5]]
    )
  end

  it 'threads a line through a corner whose value equals the level' do
    values = [
      [-1, 1, 3],
      [-2, 0, 2],
      [-3, -1, 1]
    ]
    result = contourer_class.new(values).contour(0)

    expect(result.line_strings.length).to eq(1)
    expect(coordinates(result)).to eq([
      [[0.5, 0], [1, 1], [1.5, 2]]
    ])
  end

  it 'pairs branches at an equal-level grid vertex by tangent direction' do
    values = [
      [-1, 1, -1],
      [-1, 0, -1],
      [1, 1, -1]
    ]
    result = contourer_class.new(values).contour(0)

    expect(result.line_strings.length).to eq(2)
    expect(coordinates(result)).to contain_exactly(
      [[0.5, 0], [1, 1], [1.5, 0]],
      [[0, 1.5], [1, 1], [1.5, 2]]
    )
  end

  it 'emits an equal-level grid edge once as the plateau boundary' do
    values = [
      [0, 0],
      [1, 1]
    ]
    result = contourer_class.new(values).contour(0)

    expect(coordinates(result)).to eq([[[0, 0], [1, 0]]])
  end

  it 'attaches an interior branch to an equal-level grid edge once' do
    values = [
      [0, 0, -1],
      [1, -1, -1],
      [-1, -1, -1]
    ]
    result = contourer_class.new(values).contour(0)

    expect(result.line_strings.length).to eq(2)
    expect(coordinates(result)).to contain_exactly(
      [[0.5, 0], [0.5, 1], [0, 1.5]],
      [[0, 0], [0.5, 0], [1, 0]]
    )
  end

  it 'keeps a branch attached at the correct position on a top level edge' do
    values = [
      [1, -1, -1],
      [0, 0, -1]
    ]
    result = contourer_class.new(values).contour(0)

    expect(coordinates(result)).to contain_exactly(
      [[0.5, 1], [0.5, 0]],
      [[0, 1], [0.5, 1], [1, 1]]
    )
  end

  it 'does not connect across cells containing missing values' do
    values = [
      [-1, -1, -1],
      [-1, 1, nil],
      [-1, nil, nil]
    ]
    result = contourer_class.new(values).contour(0.5)

    expect(result.line_strings.length).to eq(1)
    expect(coordinates(result)).to eq([
      [[1, 0.75], [0.75, 1]]
    ])
  end

  it 'supports a configured missing value marker' do
    values = [
      [-1, -1, -999],
      [1, 1, -999]
    ]
    result = contourer_class.new(values, missing_value: -999).contour(0)

    expect(result.line_strings.length).to eq(1)
    expect(coordinates(result)).to eq([[[1, 0.5], [0, 0.5]]])
  end

  it 'treats NaN cells as missing and honors grid coordinates' do
    values = [
      [-1, -1, -1],
      [1, 1, Float::NAN]
    ]
    result = contourer_class.new(values, 10, 20, 2, 4).contour(0)

    expect(result.line_strings.length).to eq(1)
    expect(coordinates(result)).to eq([[[12, 22.0], [10, 22.0]]])
  end
end
