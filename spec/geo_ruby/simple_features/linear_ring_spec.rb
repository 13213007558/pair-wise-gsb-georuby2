require File.expand_path(File.dirname(__FILE__) + '/../../spec_helper')

describe GeoRuby::SimpleFeatures::LinearRing do

  it 'should instantiate' do
    lr = GeoRuby::SimpleFeatures::LinearRing.new(4326)
    expect(lr).to be_instance_of(GeoRuby::SimpleFeatures::LinearRing)
  end

  describe 'Instance' do

    let(:lr) { GeoRuby::SimpleFeatures::LinearRing.from_coordinates([[10, 10], [20, 45], [45, 10], [10, 10]], 256) }

    it 'should test if contains a point' do
      expect(lr.contains_point?(GeoRuby::SimpleFeatures::Point.from_x_y(21, 21))).to be_truthy
    end

    it 'should test if not contains a point' do
      expect(lr.contains_point?(GeoRuby::SimpleFeatures::Point.from_x_y(21, 51))).to be_falsey
    end

  end

  describe 'spherical area' do
    it 'rejects adjacent antipodal points' do
      ring = described_class.from_coordinates([[0, 0], [180, 0], [0, 10], [0, 0]])

      expect { ring.spherical_area(1) }.to raise_error(ArgumentError)
    end

    it 'returns zero for a ring with fewer than three distinct consecutive vertices' do
      ring = described_class.from_coordinates([[0, 0], [0, 0], [10, 0], [0, 0]])

      expect(ring.spherical_area(1)).to be_zero
    end

    it 'treats consecutive pole points with different longitudes as duplicates' do
      ring = described_class.from_coordinates([[0, 90], [10, 90], [10, 0], [0, 90]])

      expect(ring.spherical_area(1)).to be_zero
    end
  end

end
