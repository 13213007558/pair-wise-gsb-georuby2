# -*- coding: utf-8 -*-
require File.expand_path(File.dirname(__FILE__) + '/../../spec_helper')

describe GeoRuby::SimpleFeatures::Geodesic do
  Point = GeoRuby::SimpleFeatures::Point
  Ellipsoid = GeoRuby::SimpleFeatures::Ellipsoid

  describe 'ellipsoids' do
    it 'uses WGS84 by default' do
      expect(Ellipsoid::WGS84.semi_major_axis).to eq(6_378_137.0)
      expect(Ellipsoid::WGS84.semi_minor_axis)
        .to be_within(1e-12).of(6_356_752.314245179)
    end

    it 'accepts caller-provided ellipsoid parameters' do
      ellipsoid = Ellipsoid.new(
        semi_major_axis: 6_378_137.0,
        semi_minor_axis: 6_356_752.314245
      )
      expect(ellipsoid.flattening).to be_within(1e-13).of(1.0 / 298.257223563)
    end

    it 'rejects invalid ellipsoids' do
      expect { Ellipsoid.new(semi_major_axis: 0, semi_minor_axis: 1) }
        .to raise_error(ArgumentError)
      expect { Ellipsoid.new(semi_major_axis: 1, semi_minor_axis: 2) }
        .to raise_error(ArgumentError)
      expect { Ellipsoid.new(semi_major_axis: Float::NAN, semi_minor_axis: 1) }
        .to raise_error(ArgumentError)
    end
  end

  describe 'inverse solutions' do
    it 'matches the published Flinders Peak to Buninyong Vincenty case' do
      start_point = Point.from_x_y(144.424867889, -37.9510334167)
      end_point = Point.from_x_y(143.926495528, -37.6528211389)
      result = start_point.geodesic_inverse(end_point)

      expect(result.status).to eq(:ok)
      expect(result).to be_converged
      expect(result.distance).to be_within(1e-6).of(54_972.271133)
      expect(result.initial_forward_azimuth).to be_within(2e-9).of(306.86815921)
      expect(result.final_forward_azimuth).to be_within(5e-9).of(127.17363064)
      expect(result.reverse_azimuth).to be_within(3e-9).of(307.17363064)
      expect(result.longitude_error.abs).to be < 1e-10
    end

    it 'solves a published polar meridian quadrant' do
      result = Point.from_x_y(0, 0).geodesic_inverse(Point.from_x_y(0, 90))

      expect(result.status).to eq(:ok)
      expect(result.distance).to be_within(1e-4).of(10_001_965.7293)
      expect(result.initial_forward_azimuth).to eq(0.0)
      expect(result.final_forward_azimuth).to eq(0.0)
      expect(result.reverse_azimuth).to eq(180.0)
    end

    it 'converges for the documented near-antipodal WGS84 stress case' do
      # Near-antipodal sample used by GeographicLib/Vincenty comparisons.
      start_point = Point.from_x_y(0, 0)
      end_point = Point.from_x_y(179.999999, -0.000001)
      result = start_point.geodesic_inverse(end_point)

      expect(result.status).to eq(:ok)
      expect(result).to be_converged
      expect(result.distance).to be_within(1e-3).of(20_003_931.348)
      expect(result.initial_forward_azimuth).to be_within(1e-7).of(179.999906)
      expect(result.final_forward_azimuth).to be_within(1e-7).of(180.000094)
      expect(result.iterations).to be > 20
      expect(result.distance_error.abs).to be < 1e-3
    end

    it 'normalizes crossings of the antimeridian' do
      result = Point.from_x_y(179, 0).geodesic_inverse(Point.from_x_y(-179, 0))

      expect(result.status).to eq(:ok)
      expect(result.distance).to be_within(1e-6).of(222_638.9815865)
      expect(result.initial_forward_azimuth).to eq(90.0)
      expect(result.reverse_azimuth).to eq(90.0)
    end

    it 'reports a determinate distance but nil azimuths for identical points' do
      start_point = Point.from_x_y(12, 34)
      result = start_point.geodesic_inverse(Point.from_x_y(-348, 34))

      expect(result.status).to eq(:same_point)
      expect(result.distance).to eq(0.0)
      expect(result.initial_forward_azimuth).to be_nil
      expect(result.final_forward_azimuth).to be_nil
      expect(result.message).to include('azimuth is undefined')
    end

    it 'treats alternative longitudes at a pole as the same point' do
      result = Point.from_x_y(10, 90).geodesic_inverse(Point.from_x_y(-50, 90))

      expect(result.status).to eq(:same_point)
      expect(result.distance).to eq(0.0)
      expect(result.initial_forward_azimuth).to be_nil
    end

    it 'reports ambiguous azimuth for exact antipodal points' do
      result = Point.from_x_y(10, 45).geodesic_inverse(Point.from_x_y(-170, -45))

      expect(result.status).to eq(:ambiguous_azimuth)
      expect(result).to be_ambiguous_azimuth
      expect(result.distance).to be_within(1e-4).of(20_003_931.4586)
      expect(result.initial_forward_azimuth).to be_nil
      expect(result.final_forward_azimuth).to be_nil
      expect(result.message).to include('more than one shortest geodesic')
    end

    it 'rejects non-finite coordinates and out-of-range latitude' do
      expect { Point.from_x_y(Float::NAN, 0).geodesic_inverse(Point.from_x_y(1, 1)) }
        .to raise_error(ArgumentError)
      expect { Point.from_x_y(0, 91).geodesic_inverse(Point.from_x_y(1, 1)) }
        .to raise_error(ArgumentError)
      expect { Point.from_x_y(0, 0).geodesic_inverse(Point.from_x_y(Float::INFINITY, 1)) }
        .to raise_error(ArgumentError)
    end
  end

  describe 'direct solutions and closure' do
    it 'matches the published Flinders Peak direct case' do
      start_point = Point.from_x_y(144.424867889, -37.9510334167)
      result = start_point.geodesic_direct(306.86815921, 54_972.271)

      expect(result.status).to eq(:ok)
      expect(result.endpoint.y).to be_within(1e-9).of(-37.652821139)
      expect(result.endpoint.x).to be_within(2e-9).of(143.926495528)
      expect(result.final_forward_azimuth).to be_within(5e-8).of(127.17363064)
    end

    it 'closes ordinary, polar and antimeridian geodesics' do
      cases = [
        [[144.424867889, -37.9510334167], 306.86815921, 54_972.271],
        [[0, 0], 0.0, 10_001_965.7293],
        [[0, 0], 90.0, 10_000_000.0],
        [[179.0, 0], 90.0, 222_638.98],
        [[0, 89.9999], 180.0, 20_003_931.4]
      ]

      cases.each do |coordinates, azimuth, distance|
        start_point = Point.from_x_y(*coordinates)
        direct = start_point.geodesic_direct(azimuth, distance)
        inverse = start_point.geodesic_inverse(direct.endpoint)

        expect(direct.status).to eq(:ok)
        expect(inverse.status).to eq(:ok)
        expect(inverse.distance).to be_within(1e-3).of(distance)
        expect(inverse.initial_forward_azimuth)
          .to be_within(1e-6).of(azimuth % 360.0)
      end
    end

    it 'returns the endpoint on the antimeridian in canonical longitude' do
      result = Point.from_x_y(179, 0).geodesic_direct(90, 222_638.98)

      expect(result.status).to eq(:ok)
      expect(result.endpoint.x).to be_within(1e-7).of(-179.0)
    end

    it 'returns undefined azimuths for a zero-distance direct solution' do
      start_point = Point.from_x_y(1, 2)
      result = start_point.geodesic_direct(45, 0)

      expect(result.status).to eq(:same_point)
      expect(result.endpoint).to eq(start_point)
      expect(result.initial_forward_azimuth).to be_nil
      expect(result.final_forward_azimuth).to be_nil
    end

    it 'rejects an invalid direct azimuth, distance or ellipsoid' do
      expect { Point.from_x_y(0, 0).geodesic_direct(Float::NAN, 1) }
        .to raise_error(ArgumentError)
      expect { Point.from_x_y(0, 0).geodesic_direct(0, -1) }
        .to raise_error(ArgumentError)
      expect { Point.from_x_y(0, 0).geodesic_direct(0, 1, :wgs84) }
        .to raise_error(ArgumentError)
    end
  end

  describe 'legacy planar distance' do
    it 'keeps the euclidean definition unchanged' do
      start_point = Point.from_x_y(0, 0)
      end_point = Point.from_x_y(3, 4)

      expect(start_point.euclidian_distance(end_point)).to eq(5.0)
    end
  end
end
