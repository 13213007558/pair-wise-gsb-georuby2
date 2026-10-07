require File.expand_path(File.dirname(__FILE__) + '/../../spec_helper')

describe GeoRuby::SimpleFeatures::LinearReferencing do

  # A short measured line whose M rises, then falls, then stays flat:
  #
  #   M
  #  20         C
  #  10      B     D---E---F
  #   0   A
  #       0  10  20  30  40  50   x
  #
  # M 10 occurs at vertex B, at vertex D and along the whole D-F plateau.
  let(:line) do
    GeoRuby::SimpleFeatures::LineString.from_coordinates(
      [[0, 0, 0], [10, 0, 10], [20, 0, 20], [30, 0, 10], [40, 0, 10], [50, 0, 10]],
      4326, false, true)
  end

  def m_values(geometry)
    geometry.points.map(&:m)
  end

  describe '#locate_along' do
    it 'returns a Point for a single isolated crossing' do
      result = line.locate_along(5)
      expect(result).to be_instance_of(GeoRuby::SimpleFeatures::Point)
      expect(result.x).to eql(5.0)
      expect(result.m).to eql(5.0)
    end

    it 'returns a MultiPoint when a backtrack crosses the measure twice' do
      result = line.locate_along(15)
      expect(result).to be_instance_of(GeoRuby::SimpleFeatures::MultiPoint)
      expect(result.length).to eql(2)
      expect(result.points.map(&:x)).to eql([15.0, 25.0])
      expect(result.points.map(&:m)).to eql([15.0, 15.0])
    end

    it 'distinguishes an isolated point from a continuous plateau segment' do
      result = line.locate_along(10)
      expect(result).to be_instance_of(GeoRuby::SimpleFeatures::GeometryCollection)
      point, plateau = result.geometries
      expect(point).to be_instance_of(GeoRuby::SimpleFeatures::Point)
      expect(point.x).to eq(10)
      expect(point.m).to eq(10)
      expect(plateau).to be_instance_of(GeoRuby::SimpleFeatures::LineString)
      expect(plateau.points.map(&:x)).to eq([30, 40, 50])
      expect(plateau.m_range).to eq([10, 10])
    end

    it 'reports a vertex shared by two segments only once' do
      result = line.locate_along(20)
      expect(result).to be_instance_of(GeoRuby::SimpleFeatures::Point)
      expect(result.x).to eq(20)
    end

    it 'returns nil when the measure is not on the line' do
      expect(line.locate_along(25)).to be_nil
    end

    it 'keeps the SRID on the result' do
      measured = GeoRuby::SimpleFeatures::LineString.from_coordinates(
        [[0, 0, 0], [10, 0, 10]], 3857, false, true)
      expect(measured.locate_along(5).srid).to eql(3857)
    end

    it 'interpolates Z with the same parameter as X/Y' do
      measured = GeoRuby::SimpleFeatures::LineString.from_coordinates(
        [[0, 0, 100, 0], [10, 0, 200, 10]], 4326, true, true)
      result = measured.locate_along(5)
      expect(result.z).to eql(150.0)
      expect(result.with_z).to be_truthy
      expect(result.with_m).to be_truthy
    end

    it 'raises when the geometry has no M ordinate' do
      plain = GeoRuby::SimpleFeatures::LineString.from_coordinates([[0, 0], [1, 1]])
      expect { plain.locate_along(0) }.to raise_error(ArgumentError)
    end
  end

  describe '#locate_between' do
    it 'clips with per-segment M interpolation, not line-length ratio' do
      result = line.locate_between(12, 18)
      expect(result).to be_instance_of(GeoRuby::SimpleFeatures::MultiLineString)
      expect(result.length).to eql(2)
      # rising side: M 12..18 is x 12..18; falling side: x 22..28
      expect(result[0].points.map(&:x)).to eql([12.0, 18.0])
      expect(result[1].points.map(&:x)).to eql([22.0, 28.0])
    end

    it 'checks the M range of every clipped result, segment by segment' do
      result = line.locate_between(12, 18)
      result.each do |part|
        expect(part.m_range).to eql([12.0, 18.0])
        part.points.each_cons(2) do |p1, p2|
          expect(p1.m).to be >= 12.0
          expect(p1.m).to be <= 18.0
          expect(p2.m).to be >= 12.0
          expect(p2.m).to be <= 18.0
        end
      end
    end

    it 'merges connected pieces across junctions into one line' do
      result = line.locate_between(10, 20)
      expect(result).to be_instance_of(GeoRuby::SimpleFeatures::LineString)
      expect(result.points.map(&:x)).to eq([10, 20, 30, 40, 50])
      expect(result.m_range).to eq([10, 20])
    end

    it 'keeps the original direction of travel' do
      result = line.locate_between(12, 18)
      expect(result[1].points.map(&:x)).to eql(result[1].points.map(&:x).sort)
      expect(m_values(result[1])).to eql([18.0, 12.0])
    end

    it 'accepts the interval bounds in any order' do
      expect(line.locate_between(18, 12)).to eq(line.locate_between(12, 18))
    end

    it 'equals locate_along when both bounds are the same measure' do
      expect(line.locate_between(10, 10)).to eq(line.locate_along(10))
    end

    it 'returns nil when nothing falls inside the interval' do
      expect(line.locate_between(25, 30)).to be_nil
    end

    context 'with open interval endpoints' do
      it 'drops plateau segments sitting on an open endpoint' do
        result = line.locate_between(10, 20, include_start: false)
        expect(result).to be_instance_of(GeoRuby::SimpleFeatures::LineString)
        expect(result.points.map(&:x)).to eq([10, 20, 30])
      end

      it 'drops isolated points sitting on an open endpoint' do
        closed = line.locate_between(5, 10)
        expect(closed).to be_instance_of(GeoRuby::SimpleFeatures::MultiLineString)
        expect(closed.length).to eql(2)
        opened = line.locate_between(5, 10, include_end: false)
        expect(opened).to be_instance_of(GeoRuby::SimpleFeatures::LineString)
        expect(opened.points.map(&:x)).to eq([5.0, 10])
      end
    end

    it 'keeps SRID and M flags on multi line results' do
      measured = GeoRuby::SimpleFeatures::LineString.from_coordinates(
        [[0, 0, 0], [10, 0, 10], [20, 0, 20], [30, 0, 10]], 3857, false, true)
      result = measured.locate_between(12, 18)
      expect(result.srid).to eql(3857)
      expect(result.with_m).to be_truthy
      result.each { |part| expect(part.srid).to eql(3857) }
    end
  end

  describe 'degenerate and dirty input' do
    it 'skips segments with a missing M' do
      measured = GeoRuby::SimpleFeatures::LineString.from_coordinates(
        [[0, 0, 0], [10, 0, 10], [20, 0, 20], [30, 0, 30]], 4326, false, true)
      measured[1].m = nil
      expect(measured.locate_along(5)).to be_nil
      result = measured.locate_between(0, 30)
      expect(result).to be_instance_of(GeoRuby::SimpleFeatures::LineString)
      expect(result.points.map(&:x)).to eq([20, 30])
    end

    it 'collapses duplicate vertices instead of emitting empty pieces' do
      measured = GeoRuby::SimpleFeatures::LineString.from_coordinates(
        [[0, 0, 0], [10, 0, 10], [10, 0, 10], [20, 0, 20]], 4326, false, true)
      result = measured.locate_between(0, 20)
      expect(result).to be_instance_of(GeoRuby::SimpleFeatures::LineString)
      expect(result.points.map(&:x)).to eq([0, 10, 20])
      expect(measured.locate_along(10)).to be_instance_of(GeoRuby::SimpleFeatures::Point)
    end

    it 'treats a zero-length vertex pair as a point, not a line' do
      measured = GeoRuby::SimpleFeatures::LineString.from_coordinates(
        [[0, 0, 0], [10, 0, 5], [10, 0, 5], [20, 0, 10]], 4326, false, true)
      result = measured.locate_along(5)
      expect(result).to be_instance_of(GeoRuby::SimpleFeatures::Point)
      expect(result.x).to eq(10)
    end
  end

  describe 'on a MultiLineString' do
    let(:multi) do
      parts = [
        GeoRuby::SimpleFeatures::LineString.from_coordinates(
          [[0, 0, 0], [10, 0, 10]], 4326, false, true),
        GeoRuby::SimpleFeatures::LineString.from_coordinates(
          [[10, 0, 10], [20, 0, 20]], 4326, false, true)
      ]
      multi = GeoRuby::SimpleFeatures::MultiLineString.from_line_strings(parts, 4326, false, true)
      multi.with_m = true
      multi
    end

    it 'merges parts that meet on a junction point' do
      result = multi.locate_between(0, 20)
      expect(result).to be_instance_of(GeoRuby::SimpleFeatures::LineString)
      expect(result.points.map(&:x)).to eq([0, 10, 20])
    end

    it 'returns a MultiLineString for disconnected parts' do
      parts = [
        GeoRuby::SimpleFeatures::LineString.from_coordinates(
          [[0, 0, 0], [10, 0, 10]], 4326, false, true),
        GeoRuby::SimpleFeatures::LineString.from_coordinates(
          [[100, 0, 10], [110, 0, 20]], 4326, false, true)
      ]
      multi = GeoRuby::SimpleFeatures::MultiLineString.from_line_strings(parts, 4326, false, true)
      multi.with_m = true
      result = multi.locate_between(5, 15)
      expect(result).to be_instance_of(GeoRuby::SimpleFeatures::MultiLineString)
      expect(result.length).to eql(2)
      result.each { |part| expect(part.m_range[0]).to be >= 5.0 }
    end
  end
end
