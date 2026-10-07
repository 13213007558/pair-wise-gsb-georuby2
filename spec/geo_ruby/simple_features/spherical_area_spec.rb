require File.expand_path(File.dirname(__FILE__) + '/../../spec_helper')

# Independent expectations are derived from spherical trigonometry
# (law of cosines / spherical excess), not from the implementation
# under test. No GIS software is involved.
describe GeoRuby::SimpleFeatures::SphericalArea do
  let(:radius) { 6_371_000.0 }

  def ring(coords)
    GeoRuby::SimpleFeatures::LinearRing.from_coordinates(coords)
  end

  def polygon(*rings)
    GeoRuby::SimpleFeatures::Polygon.from_coordinates(rings)
  end

  def radians(degrees)
    degrees * Math::PI / 180.0
  end

  # Angle at (lon, lat) between the meridian towards the north pole
  # and the great circle towards (lon + dlon, lat), from the spherical
  # law of cosines on the pole-vertex triangle.
  def beta(lat, dlon)
    sin_lat = Math.sin(radians(lat))
    cos_lat = Math.cos(radians(lat))
    cos_c = sin_lat**2 + cos_lat**2 * Math.cos(radians(dlon))
    cos_beta = Math.tan(radians(lat)) * Math.tan(Math.acos(cos_c) / 2)
    Math.acos(cos_beta)
  end

  # Exact area of a ring [[lon1, lat1], [lon2, lat1], [lon2, lat2],
  # [lon1, lat2]] whose edges are great-circle arcs (the top and
  # bottom edges bulge poleward, they are not parallels).
  def rectangle_ring_area(lon1, lat1, lon2, lat2, r)
    dlon = lon2 - lon1
    2 * (beta(lat1, dlon) - beta(lat2, dlon)) * r * r
  end

  describe 'LinearRing spherical area' do
    it 'computes one eighth of the sphere (equator and meridians)' do
      octant = ring([[0, 0], [90, 0], [0, 90]])
      expect(octant.spherical_area(1.0)).to be_within(1e-9).of(Math::PI / 2)
    end

    it 'computes a zone bounded by the equator and two meridians' do
      zone = ring([[0, 0], [90, 0], [90, 45], [0, 45]])
      expected = rectangle_ring_area(0, 0, 90, 45, 1.0)
      expect(zone.spherical_area(1.0)).to be_within(1e-9).of(expected)
    end

    it 'computes a lune bounded by two meridians' do
      lune = ring([[0, 0], [0, -90], [90, 0], [0, 90]])
      expect(lune.spherical_area(1.0)).to be_within(1e-9).of(Math::PI)
    end

    it 'computes a hemisphere bounded by the equator' do
      hemisphere = ring([[0, 0], [120, 0], [-120, 0]])
      expect(hemisphere.spherical_area(1.0)).to be_within(1e-9).of(2 * Math::PI)
    end

    it 'scales with the square of the caller-supplied radius' do
      zone = ring([[0, 0], [90, 0], [90, 45], [0, 45]])
      expected = rectangle_ring_area(0, 0, 90, 45, radius)
      expect(zone.spherical_area(radius)).to be_within(1e-6 * expected).of(expected)
    end

    it 'negates the signed area when the ring direction is flipped' do
      coords = [[0, 0], [90, 0], [90, 45], [0, 45]]
      expected = -rectangle_ring_area(0, 0, 90, 45, 1.0)
      expect(ring(coords.reverse).spherical_area(1.0)).to be_within(1e-9).of(expected)
    end

    it 'follows the short arc across the antimeridian' do
      zone = ring([[170, 10], [-170, 10], [-170, 20], [170, 20]])
      expected = rectangle_ring_area(170, 10, 190, 20, 1.0)
      expect(zone.spherical_area(1.0)).to be_within(1e-9).of(expected)
    end

    it 'crosses the antimeridian along the equator' do
      lune = ring([[170, 0], [-170, 0], [0, 90]])
      expect(lune.spherical_area(1.0)).to be_within(1e-9).of(radians(20))
    end

    it 'encloses the north pole' do
      cap = ring([[0, 80], [120, 80], [-120, 80]])
      expected = 6 * beta(80, 120) - Math::PI
      expect(cap.spherical_area(1.0)).to be_within(1e-9).of(expected)
    end

    it 'ignores a repeated closing point' do
      open_ring = ring([[0, 0], [90, 0], [90, 45], [0, 45]])
      closed_ring = ring([[0, 0], [90, 0], [90, 45], [0, 45], [0, 0]])
      expect(closed_ring.spherical_area(1.0))
        .to be_within(1e-12).of(open_ring.spherical_area(1.0))
    end

    it 'returns 0.0 for zero-area rings' do
      expect(ring([[0, 0], [10, 0], [20, 0]]).spherical_area(1.0)).to eq(0.0)
      expect(ring([[5, 5], [5, 5]]).spherical_area(1.0)).to eq(0.0)
    end

    it 'rejects adjacent antipodal vertices' do
      expect { ring([[0, 0], [180, 0], [0, 45]]).spherical_area(1.0) }
        .to raise_error(ArgumentError)
    end

    it 'rejects a vertex antipodal to the fan origin' do
      expect { ring([[0, 0], [30, 30], [180, 0], [30, -30]]).spherical_area(1.0) }
        .to raise_error(ArgumentError)
    end

    it 'cancels opposite lobes of a self-intersecting ring' do
      eight = ring([[0, 0], [10, 10], [20, 0], [0, 0], [10, -10], [20, 0]])
      lobe = ring([[0, 0], [10, 10], [20, 0]]).spherical_signed_area(1.0)
      expect(lobe).not_to eq(0.0)
      expect(eight.spherical_area(1.0)).to be_within(1e-9).of(0.0)
    end
  end

  describe 'regions larger than a hemisphere' do
    # ring whose left side is the area north of latitude -1 (> hemisphere)
    let(:big_ring) { ring([[0, -1], [120, -1], [-120, -1]]) }
    let(:big_area) { 6 * beta(-1, 120) - Math::PI }

    it 'reports the small complement by default (signed convention)' do
      expect(big_ring.spherical_area(1.0))
        .to be_within(1e-9).of(big_area - 4 * Math::PI)
    end

    it 'reports the large region only through an explicit rule' do
      expect(big_ring.spherical_area(1.0, rule: :left))
        .to be_within(1e-9).of(big_area)
      expect(big_ring.spherical_area(1.0, rule: :smaller))
        .to be_within(1e-9).of(4 * Math::PI - big_area)
      expect(big_ring.spherical_area(1.0, rule: :signed))
        .to be_within(1e-9).of(big_area - 4 * Math::PI)
    end

    it 'requires an explicit rule when the magnitude exceeds a hemisphere' do
      wound = ring([[0, 0], [120, 0], [-120, 0], [0, 0], [120, 0], [-120, 0]])
      expect { wound.spherical_area(1.0) }.to raise_error(ArgumentError)
      expect(wound.spherical_area(1.0, rule: :signed))
        .to be_within(1e-9).of(4 * Math::PI)
    end

    it 'rejects unknown rules' do
      expect { big_ring.spherical_area(1.0, rule: :bogus) }
        .to raise_error(ArgumentError)
    end
  end

  describe 'Polygon spherical area' do
    let(:outer) { [[0, 0], [10, 0], [10, 10], [0, 10]] }
    let(:hole_ccw) { [[2, 2], [8, 2], [8, 8], [2, 8]] }
    let(:hole_cw) { [[2, 2], [2, 8], [8, 8], [8, 2]] }

    it 'subtracts holes wound opposite to the outer ring' do
      expected = rectangle_ring_area(0, 0, 10, 10, 1.0) -
                 rectangle_ring_area(2, 2, 8, 8, 1.0)
      expect(polygon(outer, hole_cw).spherical_area(1.0))
        .to be_within(1e-9).of(expected)
    end

    it 'adds holes wound like the outer ring (signed, never absolutized)' do
      expected = rectangle_ring_area(0, 0, 10, 10, 1.0) +
                 rectangle_ring_area(2, 2, 8, 8, 1.0)
      expect(polygon(outer, hole_ccw).spherical_area(1.0))
        .to be_within(1e-9).of(expected)
    end

    it 'requires an explicit rule when the total exceeds a hemisphere' do
      wound = [[0, 0], [120, 0], [-120, 0], [0, 0], [120, 0], [-120, 0]]
      expect { polygon(wound).spherical_area(1.0) }.to raise_error(ArgumentError)
      expect(polygon(wound).spherical_area(1.0, rule: :signed))
        .to be_within(1e-9).of(4 * Math::PI)
    end
  end
end
