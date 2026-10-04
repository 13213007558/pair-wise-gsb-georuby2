require File.expand_path(File.dirname(__FILE__) + '/../../spec_helper')

describe GeoRuby::SimpleFeatures::Contour do
  describe '.from_grid' do
    it 'builds an open line from a linear grid boundary' do
      result = described_class.from_grid(
        [[0, 1], [0, 1]],
        level: 0.5,
        origin_x: 10,
        origin_y: 20,
        cell_width: 2,
        cell_height: 4,
        srid: 4326,
        curve_segments: 1
      )

      expect(result.open_lines.length).to eq(1)
      expect(result.closed_rings).to be_empty
      expect(result.open_lines.first).to be_a(GeoRuby::SimpleFeatures::LineString)
      expect(result.open_lines.first.srid).to eq(4326)
      expect(result.open_lines.first.to_coordinates).to eq([[11.0, 20.0], [11.0, 24.0]])
    end

    it 'reuses intersection identity across neighboring cells' do
      grid = [[0, 1, 0], [0, 0, 0]]
      builder = described_class.new(
        grid, 0.5, 0.0, 0.0, 1.0, 1.0,
        GeoRuby::SimpleFeatures::DEFAULT_SRID, nil, 1
      )
      result = builder.build

      expect(result.open_lines.length).to eq(1)
      line = result.open_lines.first
      expect(line.to_coordinates).to eq([[0.5, 0.0], [1.0, 0.5], [1.5, 0.0]])
      expect(line.points[1]).to be_a(GeoRuby::SimpleFeatures::Point)
      shared_point = builder.instance_variable_get(:@points).values.find do |point|
        point.x == 1.0 && point.y == 0.5
      end
      expect(line.points[1]).to equal(shared_point)
    end

    it 'samples bilinear curves inside cells' do
      result = described_class.from_grid(
        [[0, 1, 0], [0, 0, 0]],
        level: 0.5,
        curve_segments: 4
      )

      expect(result.open_lines.length).to eq(1)
      line = result.open_lines.first
      expect(line.length).to eq(9)
      expect(line.first.to_coordinates).to eq([0.5, 0.0])
      expect(line.last.to_coordinates).to eq([1.5, 0.0])
      expect(line.points[1..-2].map(&:to_coordinates)).to all(
        satisfy { |coordinate| coordinate[0].between?(0.5, 1.5) }
      )
    end

    it 'builds a closed linear ring around a high region' do
      result = described_class.from_grid(
        [[0, 1, 0], [1, 2, 1], [0, 1, 0]],
        level: 1,
        curve_segments: 1
      )

      expect(result.open_lines).to be_empty
      expect(result.closed_rings.length).to eq(1)
      ring = result.closed_rings.first
      expect(ring).to be_a(GeoRuby::SimpleFeatures::LinearRing)
      expect(ring).to be_closed
      expect(ring.to_coordinates).to eq(
        [[1.0, 0.0], [0.0, 1.0], [1.0, 2.0], [2.0, 1.0], [1.0, 0.0]]
      )
    end

    it 'wraps closed rings in polygons on demand' do
      result = described_class.from_grid(
        [[0, 1, 0], [1, 2, 1], [0, 1, 0]],
        level: 1,
        curve_segments: 1
      )

      expect(result.polygons.length).to eq(1)
      expect(result.polygons.first).to be_a(GeoRuby::SimpleFeatures::Polygon)
    end

    it 'resolves a saddle below the level with disconnected diagonals' do
      result = described_class.from_grid(
        [[-1, 1], [1, -0.25]],
        level: 0,
        curve_segments: 1
      )

      expect(result.closed_rings).to be_empty
      expect(result.open_lines.map(&:to_coordinates)).to match_array([
        [[0.5, 0.0], [0.0, 0.5]],
        [[1.0, 0.8], [0.8, 1.0]]
      ])
    end

    it 'changes saddle topology when the level moves above the center' do
      result = described_class.from_grid(
        [[-1, 1], [1, 0.25]],
        level: 0,
        curve_segments: 1
      )

      expect(result.closed_rings).to be_empty
      expect(result.open_lines.map(&:to_coordinates)).to eq([
        [[0.5, 0.0], [0.0, 0.5]]
      ])
    end

    it 'emits two crossing open lines when the level passes exactly through a saddle' do
      result = described_class.from_grid(
        [[-1, 1], [1, -1]],
        level: 0
      )

      expect(result.closed_rings).to be_empty
      expect(result.open_lines.map(&:to_coordinates)).to match_array([
        [[0.0, 0.5], [0.5, 0.5], [1.0, 0.5]],
        [[0.5, 0.0], [0.5, 0.5], [0.5, 1.0]]
      ])
      center_points = result.open_lines.map { |line| line.points[1] }
      expect(center_points.first).to eq(center_points.last)
    end

    it 'uses a grid corner equal to the level as a shared crossing point' do
      result = described_class.from_grid(
        [[0, 1], [-1, 1]],
        level: 0,
        curve_segments: 1
      )

      expect(result.open_lines.map(&:to_coordinates)).to eq([
        [[0.0, 0.0], [0.5, 1.0]]
      ])
    end

    it 'does not emit geometry for a cell entirely at the level' do
      result = described_class.from_grid(
        [[0, 0], [0, 0]],
        level: 0
      )

      expect(result.open_lines).to be_empty
      expect(result.closed_rings).to be_empty
    end

    it 'skips cells that contain nil values and ends contours at their boundary' do
      result = described_class.from_grid(
        [[0, 1, nil], [0, 1, nil]],
        level: 0.5,
        curve_segments: 1
      )

      expect(result.open_lines.map(&:to_coordinates)).to eq([
        [[0.5, 0.0], [0.5, 1.0]]
      ])
    end

    it 'skips cells containing the configured nodata value' do
      result = described_class.from_grid(
        [[0, 1, -9999], [0, 1, -9999]],
        level: 0.5,
        nodata: -9999,
        curve_segments: 1
      )

      expect(result.open_lines.map(&:to_coordinates)).to eq([
        [[0.5, 0.0], [0.5, 1.0]]
      ])
    end

    it 'treats NaN as missing data' do
      result = described_class.from_grid(
        [[0, 1, Float::NAN], [0, 1, Float::NAN]],
        level: 0.5,
        curve_segments: 1
      )

      expect(result.open_lines.map(&:to_coordinates)).to eq([
        [[0.5, 0.0], [0.5, 1.0]]
      ])
    end

    it 'returns no geometry outside the value range' do
      result = described_class.from_grid(
        [[0, 1], [0, 1]],
        level: 2
      )

      expect(result.open_lines).to be_empty
      expect(result.closed_rings).to be_empty
    end
  end
end
