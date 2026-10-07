# -*- coding: utf-8 -*-
require File.expand_path(File.dirname(__FILE__) + '/../../spec_helper')

describe GeoRuby::SimpleFeatures::Geodesic do
  let(:wgs84) { GeoRuby::SimpleFeatures::Ellipsoid::WGS84 }

  def point(lon, lat)
    GeoRuby::SimpleFeatures::Point.from_x_y(lon, lat)
  end

  # Independent reference: integrates the exact geodesic ODEs on the
  # ellipsoid (meridian radius of curvature, prime-vertical radius and
  # Clairaut's relation) with a classical RK4 scheme.  Shares no code
  # with the Fourier-series implementation under test.
  def rk4_direct(lat1, lon1, azi1_deg, s12, steps = 3000)
    a = 6_378_137.0
    f = 1.0 / 298.257223563
    e2 = f * (2 - f)
    deg = Math::PI / 180.0
    phi = lat1 * deg
    lam = lon1 * deg
    alpha = azi1_deg * deg
    ds = s12 / steps
    deriv = lambda do |ph, _lm, al|
      s = Math.sin(ph)
      w = Math.sqrt(1 - e2 * s * s)
      rho = a * (1 - e2) / (w**3)
      rn = a / w
      beta = Math.atan((1 - f) * Math.tan(ph))
      cb = Math.cos(beta)
      cp = Math.cos(ph)
      [Math.cos(al) / rho,
       Math.sin(al) / (rn * cp),
       Math.tan(beta) * Math.sin(al) * (1 - f) * cb * cb / (cp * cp * rho)]
    end
    steps.times do
      k1 = deriv.call(phi, lam, alpha)
      k2 = deriv.call(phi + ds * k1[0] / 2, lam + ds * k1[1] / 2, alpha + ds * k1[2] / 2)
      k3 = deriv.call(phi + ds * k2[0] / 2, lam + ds * k2[1] / 2, alpha + ds * k2[2] / 2)
      k4 = deriv.call(phi + ds * k3[0], lam + ds * k3[1], alpha + ds * k3[2])
      phi += ds * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]) / 6
      lam += ds * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]) / 6
      alpha += ds * (k1[2] + 2 * k2[2] + 2 * k3[2] + k4[2]) / 6
    end
    [phi / deg, ((lam / deg + 180) % 360) - 180]
  end

  # Meridian arc from the equator to geodetic latitude +lat_deg+, by
  # Simpson quadrature of the meridian radius of curvature.  Independent
  # of the Fourier-series code under test.
  def meridian_arc(lat_deg)
    a = 6_378_137.0
    e2 = (1.0 / 298.257223563) * (2 - 1.0 / 298.257223563)
    n = 2000
    h = lat_deg * Math::PI / 180.0 / n
    rho = lambda do |phi|
      s = Math.sin(phi)
      a * (1 - e2) / (1 - e2 * s * s)**1.5
    end
    sum = rho.call(0.0) + rho.call(lat_deg * Math::PI / 180.0)
    (1...n).each { |i| sum += (i.odd? ? 4 : 2) * rho.call(i * h) }
    sum * h / 3
  end

  def load_samples
    path = File.expand_path(File.dirname(__FILE__) +
                            '/../../../data/geodesic/geodetic_samples.dat')
    File.readlines(path).reject { |l| l.strip.empty? || l.start_with?('#') }
        .map { |l| l.split.map(&:to_f) }
  end

  describe 'published WGS84 reference constants' do
    it 'computes the meridian quadrant (equator to pole)' do
      r = point(0, 0).geodesic_inverse(point(0, 90))
      expect(r.distance).to be_within(1e-4).of(10_001_965.7293)
      expect(r.azimuth1).to eq(0.0)
      expect(r).to be_converged
    end

    it 'computes the equatorial quadrant' do
      r = point(0, 0).geodesic_inverse(point(90, 0))
      expect(r.distance).to be_within(1e-4).of(10_018_754.1714)
      expect(r.azimuth1).to be_within(1e-12).of(90.0)
      expect(r.azimuth2).to be_within(1e-12).of(90.0)
    end

    it 'computes the half meridian (equator to equator over the pole)' do
      r = point(0, 0).geodesic_inverse(point(180, 0))
      expect(r.distance).to be_within(1e-4).of(20_003_931.4586)
    end
  end

  describe 'geodetic test samples (ordinary, polar, near-antipodal)' do
    it 'matches the reference distances and closes the round trip' do
      expect(load_samples.size).to be >= 10
      load_samples.each do |lat1, lon1, lat2, lon2, s12|
        p1 = point(lon1, lat1)
        p2 = point(lon2, lat2)
        r = p1.geodesic_inverse(p2)
        expect(r).to be_converged
        expect(r.distance).to be_within(2e-3).of(s12)
        expect(r.error_estimate).to be < 1e-4
        # direct/inverse closure: flying azimuth1 for distance lands on p2
        d = p1.geodesic_direct(r.azimuth1, r.distance)
        expect(d.point.y).to be_within(1e-9).of(lat2)
        dlon = (d.point.x - lon2 + 540) % 360 - 180
        expect(dlon.abs).to be < 1e-9
        # azimuth2 is consistent with the reverse problem
        back = p2.geodesic_inverse(p1)
        expect(back.distance).to be_within(1e-6).of(r.distance)
        dazi = (back.azimuth1 - (r.azimuth2 + 180) % 360).abs
        dazi = 360 - dazi if dazi > 180
        expect(dazi).to be < 1e-7
      end
    end
  end

  describe 'independent numerical cross-check' do
    it 'agrees with RK4 integration of the exact geodesic ODEs' do
      [
        [40.6, -73.8, 53.4, 5_500_000.0],   # ordinary
        [80.0, 10.0, 170.0, 4_000_000.0],   # polar start
        [30.0, 0.0, 29.9, 19_900_000.0],    # near-antipodal
        [10.0, 170.0, 80.0, 3_000_000.0]    # dateline crossing
      ].each do |lat1, lon1, azi, s12|
        d = point(lon1, lat1).geodesic_direct(azi, s12)
        rlat, rlon = rk4_direct(lat1, lon1, azi, s12)
        expect(d.point.y).to be_within(1e-6).of(rlat)
        dlon = (d.point.x - rlon + 540) % 360 - 180
        expect(dlon.abs).to be < 1e-6
      end
    end
  end

  describe 'near-antipodal robustness' do
    it 'converges where plain Vincenty-style iteration may not' do
      # This nearly antipodal pair is beyond the reach of Vincenty's
      # iteration (Point#ellipsoidal_distance returns NaN when it does
      # not converge); the geodesic inverse must still return a
      # converged ellipsoidal result - never a spherical approximation.
      p1 = point(-11.608970491806957, 27.042118541426206)
      p2 = point(168.38179539728412, -27.33726945112841)
      vincenty = p1.ellipsoidal_distance(p2)
      r = p1.geodesic_inverse(p2)
      expect(r).to be_converged
      expect(r.distance).to be_within(1e-3).of(19_971_221.9418)
      expect(r.azimuth1).to be_within(1e-6).of(179.390425)
      expect(r.azimuth2).to be_within(1e-6).of(0.611182)
      # when Vincenty does converge, both ellipsoidal results agree
      expect(r.distance).to be_within(1e-3).of(vincenty) if vincenty.finite?
    end

    it 'finds the shortest geodesic among several near-antipodal ones' do
      p1 = point(-11.608970491806957, 27.042118541426206)
      p2 = point(168.38179539728412, -27.33726945112841)
      r = p1.geodesic_inverse(p2)
     expect(r.distance).to be < 20_003_931.4586
   end
 end
  describe 'special configurations' do
    it 'handles coincident points: distance certain, azimuth undefined' do
      r = point(30, 45).geodesic_inverse(point(30, 45))
      expect(r).to be_converged
      expect(r.distance).to eq(0.0)
      expect(r.azimuth1).to be_nil
      expect(r.azimuth2).to be_nil
      expect(r).not_to be_azimuths_unique
    end

    it 'treats two references to the same pole as coincident' do
      r = point(0, 90).geodesic_inverse(point(120, 90))
      expect(r.distance).to eq(0.0)
      expect(r).not_to be_azimuths_unique
    end

    it 'handles pole to pole: distance certain, azimuth not unique' do
      r = point(0, 90).geodesic_inverse(point(0, -90))
      expect(r).to be_converged
      expect(r.distance).to be_within(1e-4).of(20_003_931.4586)
      expect(r).not_to be_azimuths_unique
    end

    it 'handles exactly antipodal points: distance certain, azimuth not unique' do
      r = point(0, 30).geodesic_inverse(point(180, -30))
      expect(r).to be_converged
      expect(r.distance).to be_within(1e-4).of(20_003_931.4586)
      expect(r).not_to be_azimuths_unique
    end

    it 'picks the shorter polar path on opposite meridians' do
      # (10, 0) -> (20, 180): shorter over the north pole; the expected
      # distance is two quadrant arcs minus the two meridian arcs,
      # computed here by independent numerical quadrature.
      r = point(0, 10).geodesic_inverse(point(180, 20))
      expect(r).to be_converged
      expect(r.azimuth1).to eq(0.0)
      expect(r.azimuth2).to eq(180.0)
      expect(r).to be_azimuths_unique
      expected = 2 * 10_001_965.7293 - meridian_arc(10) - meridian_arc(20)
      expect(r.distance).to be_within(1e-3).of(expected)
    end

    it 'handles points on the same meridian' do
      r = point(30, 10).geodesic_inverse(point(30, 55))
      expect(r.azimuth1).to eq(0.0)
      expect(r.azimuth2).to eq(0.0)
      r = point(30, 55).geodesic_inverse(point(30, 10))
      expect(r.azimuth1).to eq(180.0)
    end

    it 'normalizes longitude differences across the dateline' do
      r1 = point(170, 10).geodesic_inverse(point(-170, 10))
      r2 = point(170, 10).geodesic_inverse(point(190, 10))
      expect(r1.distance).to be_within(1e-9).of(r2.distance)
      expect(r1.azimuth1).to be_within(1e-9).of(r2.azimuth1)
    end


    it 'wraps the destination longitude into [-180, 180)' do
      d = point(170, 10).geodesic_direct(90, 3_000_000)
      expect(d.point.x).to be >= -180.0
      expect(d.point.x).to be < 180.0
      # closure: the inverse from the destination comes back with 3000 km
      back = d.point.geodesic_inverse(point(170, 10))
      expect(back.distance).to be_within(1e-6).of(3_000_000)
    end

    it 'departs from a pole along the conventional meridian' do
      # from the north pole, azimuth azi leaves along lon1 + 180 - azi
      d = point(0, 90).geodesic_direct(90, 1_000_000)
      expect(d.point.x).to be_within(1e-9).of(90.0)
      expect(point(0, 90).geodesic_inverse(d.point).distance)
        .to be_within(1e-6).of(1_000_000)
      # from the south pole, azimuth azi leaves along lon1 + azi
      d = point(10, -90).geodesic_direct(45, 1_000_000)
      expect(d.point.x).to be_within(1e-9).of(55.0)
      expect(point(10, -90).geodesic_inverse(d.point).distance)
        .to be_within(1e-6).of(1_000_000)
    end

    it 'supports negative distances (backwards along the geodesic)' do
      d1 = point(0, 0).geodesic_direct(90, 1_000_000)
      d2 = point(0, 0).geodesic_direct(270, -1_000_000)
      expect(d2.point.y).to be_within(1e-12).of(d1.point.y)
      expect(d2.point.x).to be_within(1e-12).of(d1.point.x)
    end

    it 'returns the origin for a zero distance' do
      d = point(30, 45).geodesic_direct(77, 0)
      expect(d.point.y).to eq(45.0)
      expect(d.point.x).to eq(30.0)
      expect(d.azimuth).to be_within(1e-12).of(77.0)
    end

    it 'accepts a custom ellipsoid as object, array or hash' do
      p1 = point(0, 10)
      p2 = point(50, 20)
      ells = [GeoRuby::SimpleFeatures::Ellipsoid.new(6_378_145.0, 1.0 / 293.0),
              [6_378_145.0, 1.0 / 293.0],
              { a: 6_378_145.0, f: 1.0 / 293.0 }]
      distances = ells.map { |e| p1.geodesic_inverse(p2, ellipsoid: e).distance }
      expect(distances[1]).to be_within(1e-9).of(distances[0])
      expect(distances[2]).to be_within(1e-9).of(distances[0])
    end

    it 'supports a spherical ellipsoid (f = 0)' do
      r = point(0, 0).geodesic_inverse(point(90, 0),
                                       ellipsoid: [6_371_000.0, 0.0])
     expect(r.distance).to be_within(1e-6).of(6_371_000.0 * Math::PI / 2)
   end
 end
  describe 'input validation' do
    it 'rejects non-finite coordinates' do
      expect { point(0, Float::NAN).geodesic_inverse(point(0, 0)) }
        .to raise_error(ArgumentError)
      expect { point(Float::INFINITY, 0).geodesic_inverse(point(0, 0)) }
        .to raise_error(ArgumentError)
      expect { point(0, 0).geodesic_inverse(point(0, Float::NAN)) }
        .to raise_error(ArgumentError)
    end

    it 'rejects latitudes outside [-90, 90]' do
      expect { point(0, 91).geodesic_inverse(point(0, 0)) }
        .to raise_error(ArgumentError)
      expect { point(0, -90.5).geodesic_direct(0, 100) }
        .to raise_error(ArgumentError)
    end

    it 'rejects non-finite azimuth or distance' do
      expect { point(0, 0).geodesic_direct(Float::NAN, 100) }
        .to raise_error(ArgumentError)
      expect { point(0, 0).geodesic_direct(45, Float::INFINITY) }
        .to raise_error(ArgumentError)
    end

    it 'rejects illegal ellipsoids' do
      expect { GeoRuby::SimpleFeatures::Ellipsoid.new(0, 0.01) }
        .to raise_error(ArgumentError)
      expect { GeoRuby::SimpleFeatures::Ellipsoid.new(-6_378_137.0, 0.01) }
        .to raise_error(ArgumentError)
      expect { GeoRuby::SimpleFeatures::Ellipsoid.new(6_378_137.0, 1.0) }
        .to raise_error(ArgumentError)
      expect { GeoRuby::SimpleFeatures::Ellipsoid.new(6_378_137.0, -1.5) }
        .to raise_error(ArgumentError)
      expect { GeoRuby::SimpleFeatures::Ellipsoid.new(6_378_137.0, Float::NAN) }
        .to raise_error(ArgumentError)
      expect { GeoRuby::SimpleFeatures::Ellipsoid.new(Float::INFINITY, 0.01) }
        .to raise_error(ArgumentError)
    end

    it 'normalizes out-of-range longitudes and azimuths instead of failing' do
      r = point(370, 10).geodesic_inverse(point(20, 20))
      expect(r).to be_converged
      d = point(0, 0).geodesic_direct(450, 1_000_000)
      expect(d.azimuth).to be >= 0.0
      expect(d.azimuth).to be < 360.0
    end
  end

  describe 'pre-existing distance methods (definitions unchanged)' do
    it 'keeps euclidian_distance planar' do
      expect(point(2, 1).euclidian_distance(point(5, 3))).to eq(Math.sqrt(13))
    end
  end
end
