require File.expand_path(File.dirname(__FILE__) + '/../../spec_helper')

module LineSpecHelper
  def mock_points(num)
    #    @point = mock(GeoRuby::SimpleFeatures::Point, :x => 1.0, :y => 2.0)
    Array.new(num) { |i| mock_point(i, i) }
  end

  def mock_point(x = 1, y = 2)
    double(GeoRuby::SimpleFeatures::Point, x: x, y: y, text_representation: "#{x} #{y}", to_coordinates: [x, y])
  end
end

describe GeoRuby::SimpleFeatures::LineString do

  include LineSpecHelper

  describe 'Instance Methods' do

    let(:line) { GeoRuby::SimpleFeatures::LineString.from_points([double(GeoRuby::SimpleFeatures::Point)]) }

    it 'should instantiate' do
      violated unless line
    end

    it 'should have binary_geometry_type 2' do
      expect(line.binary_geometry_type).to eql(2)
    end

    it 'should have text_geometry_type' do
      expect(line.text_geometry_type).to eql('LINESTRING')
    end

    it 'should have a points array' do
      expect(line.points).to be_instance_of(Array)
    end

  end

  describe 'Valid GeoRuby::SimpleFeatures::LineString' do

    let(:line) { GeoRuby::SimpleFeatures::LineString.from_points([GeoRuby::SimpleFeatures::Point.xy(1, 1), GeoRuby::SimpleFeatures::Point.xy(2, 2), GeoRuby::SimpleFeatures::Point.xy(3, 3)]) }

    it 'should check orientation' do
      expect(line).not_to be_clockwise
    end

    it 'should check orientation' do
      l = GeoRuby::SimpleFeatures::LineString.from_points([GeoRuby::SimpleFeatures::Point.from_x_y(20, 20), GeoRuby::SimpleFeatures::Point.from_x_y(10, 10), GeoRuby::SimpleFeatures::Point.from_x_y(-10, 10)])
      expect(l).to be_clockwise
    end

  end

  describe 'tu converted' do

    it 'should concat points' do
      line_string = GeoRuby::SimpleFeatures::LineString.new
      line_string.concat([GeoRuby::SimpleFeatures::Point.from_x_y(12.4, 45.3), GeoRuby::SimpleFeatures::Point.from_x_y(45.4, 41.6)])

      expect(line_string.length).to eql(2)
      expect(line_string[0]).to eq(GeoRuby::SimpleFeatures::Point.from_x_y(12.4, 45.3))

      point = GeoRuby::SimpleFeatures::Point.from_x_y(123, 45.8777)
      line_string[0] = point
      expect(line_string[0]).to eq(point)

      points = [GeoRuby::SimpleFeatures::Point.from_x_y(123, 45.8777), GeoRuby::SimpleFeatures::Point.from_x_y(45.4, 41.6)]
      line_string.each_index { |i| expect(line_string[i]).to eq(points[i]) }

      point = GeoRuby::SimpleFeatures::Point.from_x_y(22.4, 13.56)
      line_string << point
      expect(line_string.length).to eql(3)
      expect(line_string[2]).to eql(point)
    end

    it 'should create' do
      line_string = GeoRuby::SimpleFeatures::LineString.from_points([GeoRuby::SimpleFeatures::Point.from_x_y(12.4, -45.3), GeoRuby::SimpleFeatures::Point.from_x_y(45.4, 41.6)], 123)
      expect(line_string.class).to eql(GeoRuby::SimpleFeatures::LineString)
      expect(line_string.length).to eql(2)
      expect(line_string[0]).to eq(GeoRuby::SimpleFeatures::Point.from_x_y(12.4, -45.3))
      expect(line_string[1]).to eq(GeoRuby::SimpleFeatures::Point.from_x_y(45.4, 41.6))

      line_string = GeoRuby::SimpleFeatures::LineString.from_coordinates([[12.4, -45.3], [45.4, 41.6], [4.456, 1.0698]], 123)
      expect(line_string.class).to eql(GeoRuby::SimpleFeatures::LineString)
      expect(line_string.length).to eql(3)
      expect(line_string[0]).to eq(GeoRuby::SimpleFeatures::Point.from_x_y(12.4, -45.3))
      expect(line_string[1]).to eq(GeoRuby::SimpleFeatures::Point.from_x_y(45.4, 41.6))

      line_string = GeoRuby::SimpleFeatures::LineString.from_coordinates([[12.4, -45.3, 123], [45.4, 41.6, 333], [4.456, 1.0698, 987]], 123, true)
      expect(line_string.class).to eql(GeoRuby::SimpleFeatures::LineString)
      expect(line_string.length).to eql(3)
      expect(line_string[0]).to eq(GeoRuby::SimpleFeatures::Point.from_x_y_z(12.4, -45.3, 123, 123))

      line_string = GeoRuby::SimpleFeatures::LineString.from_coordinates([[12.4, -45.3, 123], [45.4, 41.6, 333], [4.456, 1.0698, 987]], 123, true)
      expect(line_string.class).to eql(GeoRuby::SimpleFeatures::LineString)
      expect(line_string.length).to eql(3)
      expect(line_string[0]).to eq(GeoRuby::SimpleFeatures::Point.from_x_y_z(12.4, -45.3, 123, 123))
    end

    it 'should bbox it' do
      bbox = GeoRuby::SimpleFeatures::LineString.from_coordinates([[12.4, -45.3, 123], [45.4, 41.6, 333], [4.456, 1.0698, 987]], 123, true).bounding_box
      expect(bbox.length).to eql(2)
      expect(bbox[0]).to eq(GeoRuby::SimpleFeatures::Point.from_x_y_z(4.456, -45.3, 123))
      expect(bbox[1]).to eq(GeoRuby::SimpleFeatures::Point.from_x_y_z(45.4, 41.6, 987))
    end

    it 'should test_line_string_equal' do
      line_string1 = GeoRuby::SimpleFeatures::LineString.from_coordinates([[12.4, -45.3], [45.4, 41.6], [4.456, 1.0698]], 123)
      line_string2 = GeoRuby::SimpleFeatures::LineString.from_coordinates([[12.4, -45.3], [45.4, 41.6]], 123)
      point = GeoRuby::SimpleFeatures::Point.from_x_y(12.4, -45.3, 123)

      expect(line_string1).to eq(GeoRuby::SimpleFeatures::LineString.from_coordinates([[12.4, -45.3], [45.4, 41.6], [4.456, 1.0698]], 123))
      expect(line_string1).not_to eq(line_string2)
      expect(line_string1).not_to eq(point)
    end

    it 'should test_line_string_binary' do
      line_string = GeoRuby::SimpleFeatures::LineString.from_coordinates([[12.4, -45.3], [45.4, 41.6]], 256)
      expect(line_string.as_hex_ewkb).to eql('01020000200001000002000000CDCCCCCCCCCC28406666666666A646C03333333333B34640CDCCCCCCCCCC4440')

      line_string = GeoRuby::SimpleFeatures::LineString.from_coordinates([[12.4, -45.3, 35.3], [45.4, 41.6, 12.3]], 256, true)
      expect(line_string.as_hex_ewkb).to eql('01020000A00001000002000000CDCCCCCCCCCC28406666666666A646C06666666666A641403333333333B34640CDCCCCCCCCCC44409A99999999992840')

      line_string = GeoRuby::SimpleFeatures::LineString.from_coordinates([[12.4, -45.3, 35.3, 45.1], [45.4, 41.6, 12.3, 40.23]], 256, true, true)
      expect(line_string.as_hex_ewkb).to eql('01020000E00001000002000000CDCCCCCCCCCC28406666666666A646C06666666666A64140CDCCCCCCCC8C46403333333333B34640CDCCCCCCCCCC44409A999999999928403D0AD7A3701D4440')
    end

    it 'test_line_string_text' do
      line_string = GeoRuby::SimpleFeatures::LineString.from_coordinates([[12.4, -45.3], [45.4, 41.6]], 256)
      expect(line_string.as_ewkt).to eql('SRID=256;LINESTRING(12.4 -45.3,45.4 41.6)')

      line_string = GeoRuby::SimpleFeatures::LineString.from_coordinates([[12.4, -45.3, 35.3], [45.4, 41.6, 12.3]], 256, true)
      expect(line_string.as_ewkt).to eql('SRID=256;LINESTRING(12.4 -45.3 35.3,45.4 41.6 12.3)')

      line_string = GeoRuby::SimpleFeatures::LineString.from_coordinates([[12.4, -45.3, 35.3], [45.4, 41.6, 12.3]], 256, false, true)
      expect(line_string.as_ewkt).to eql('SRID=256;LINESTRINGM(12.4 -45.3 35.3,45.4 41.6 12.3)')

      line_string = GeoRuby::SimpleFeatures::LineString.from_coordinates([[12.4, -45.3, 35.3, 25.2], [45.4, 41.6, 12.3, 13.75]], 256, true, true)
      expect(line_string.as_ewkt).to eql('SRID=256;LINESTRING(12.4 -45.3 35.3 25.2,45.4 41.6 12.3 13.75)')
    end

    it 'should test_linear_ring_creation' do
      # testing just the constructor helpers since the rest is the same as for line_string
      linear_ring = GeoRuby::SimpleFeatures::LinearRing.from_coordinates([[12.4, -45.3], [45.4, 41.6], [4.456, 1.0698], [12.4, -45.3]], 345)
      expect(linear_ring.class).to eql(GeoRuby::SimpleFeatures::LinearRing)
      expect(linear_ring.length).to eql(4)
      expect(linear_ring).to be_closed
      expect(linear_ring[1]).to eq(GeoRuby::SimpleFeatures::Point.from_x_y(45.4, 41.6, 345))
    end
  end

  describe '> Coordinates' do

    before(:each) do
      expect(GeoRuby::SimpleFeatures::Point).to receive(:from_coordinates)
        .exactly(4).with(anything, 4326, false, false).and_return(mock_point)
      @line = GeoRuby::SimpleFeatures::LineString.from_coordinates([1.2, 2.5, 2.2, 4.5])
    end

    it 'should instantiate from coordinates' do
      expect(@line.points.length).to eql(4)
    end

  end

  describe '> Instantiated' do

    let (:line) { GeoRuby::SimpleFeatures::LineString.from_points(mock_points(7)) }

    it 'should be closed if the last point equals the first' do
      line.push(line.first)
      expect(line).to be_closed
      expect(line.length).to eql(8)
    end

    it 'should print the text representation' do
      expect(line.text_representation).to eql('0 0,1 1,2 2,3 3,4 4,5 5,6 6')
    end

    it 'should print the georss_simple_representation' do
      expect(line.georss_simple_representation(geom_attr: nil))
        .to eql("<georss:line>0 0 1 1 2 2 3 3 4 4 5 5 6 6</georss:line>\n")
    end

    it 'should map the georss_poslist' do
      expect(line.georss_poslist).to eql('0 0 1 1 2 2 3 3 4 4 5 5 6 6')
    end

    it 'should print the kml_representation' do
      expect(line.kml_representation)
        .to eql("<LineString>\n<coordinates>0,0 1,1 2,2 3,3 4,4 5,5 6,6</coordinates>\n</LineString>\n")
    end

    it 'should print the kml_poslist without reverse or z' do
      expect(line.kml_poslist({})).to eql('0,0 1,1 2,2 3,3 4,4 5,5 6,6')
    end

    it 'should print the kml_poslist reverse' do
      expect(line.kml_poslist(reverse: true)).to eql('6,6 5,5 4,4 3,3 2,2 1,1 0,0')
    end

    it 'should print out nicely as json/geojson' do
      expect(line.as_json).to eql(type: 'LineString', coordinates: [[0, 0], [1, 1], [2, 2], [3, 3], [4, 4], [5, 5], [6, 6]])
    end

    it 'should print out nicely to json/geojson' do
      expect(line.to_json).to eql("{\"type\":\"LineString\",\"coordinates\":[[0,0],[1,1],[2,2],[3,3],[4,4],[5,5],[6,6]]}")
    end

  end

  describe '> Distances...' do
    before(:each) do
      @p1 = double(GeoRuby::SimpleFeatures::Point)
      @p2 = double(GeoRuby::SimpleFeatures::Point)
      @p3 = double(GeoRuby::SimpleFeatures::Point)
      @line = GeoRuby::SimpleFeatures::LineString.from_points([@p1, @p2, @p3])
    end

    it 'should print the length with haversine' do
      expect(@p1).to receive(:spherical_distance).with(@p2).and_return(10)
      expect(@p2).to receive(:spherical_distance).with(@p3).and_return(10)
      expect(@line.spherical_distance).to eql(20)
    end

    it 'should print lenght as euclidian' do
      expect(@p1).to receive(:euclidian_distance).with(@p2).and_return(10)
      expect(@p2).to receive(:euclidian_distance).with(@p3).and_return(10)
      expect(@line.euclidian_distance).to eql(20)
    end
  end

  describe 'Continuous Fréchet distance' do
    let(:dense_line) do
      GeoRuby::SimpleFeatures::LineString.from_coordinates([
        [0, 0], [1000, 0], [4000, 0], [8000, 0], [10000, 0]
      ])
    end

    let(:shifted_samples_line) do
      GeoRuby::SimpleFeatures::LineString.from_coordinates([
        [0, 0], [3000, 0], [6000, 0], [7000, 0], [10000, 0]
      ])
    end

    def verify_monotone_evidence(decision)
      expect(decision.within?).to be(true)
      expect(decision.evidence).not_to be_empty
      expect(decision.evidence.first[:u]).to eq(0.0)
      expect(decision.evidence.first[:v]).to eq(0.0)
      expect(decision.evidence.last[:u]).to eq(1.0)
      expect(decision.evidence.last[:v]).to eq(1.0)

      decision.evidence.each_cons(2) do |left_step, right_step|
        expect(right_step[:u]).to be >= left_step[:u]
        expect(right_step[:v]).to be >= left_step[:v]
      end

      slack = decision.threshold + decision.numeric_tolerance
      decision.evidence.each do |step|
        expect(step[:distance]).to be <= slack
        expect(step[:first_point].length).to eq(2)
        expect(step[:second_point].length).to eq(2)
      end

      decision.evidence.each_cons(2) do |left_step, right_step|
        first_span = right_step[:u] - left_step[:u]
        second_span = right_step[:v] - left_step[:v]
        largest_span = [first_span, second_span].max
        next if largest_span.zero?

        4.times do |index|
          fraction = index / 4.0
          point = left_step[:first_point].each_with_index.map do |coordinate, axis|
            coordinate + (right_step[:first_point][axis] - coordinate) * fraction
          end
          other = left_step[:second_point].each_with_index.map do |coordinate, axis|
            coordinate + (right_step[:second_point][axis] - coordinate) * fraction
          end
          distance = Math.sqrt((point[0] - other[0])**2 + (point[1] - other[1])**2)
          expect(distance).to be <= slack + 1e-10
        end
      end
    end

    it 'measures continuous movement instead of corresponding vertex distances' do
      vertex_distances = dense_line.points.zip(shifted_samples_line.points).map do |point, other|
        point.euclidian_distance(other)
      end
      expect(vertex_distances.max).to eq(2000.0)

      result = dense_line.continuous_frechet_distance(shifted_samples_line)
      decision = dense_line.continuous_frechet_threshold_decision(
        shifted_samples_line, 1e-6
      )

      expect(result.distance).to be_within(result.error_bound).of(0.0)
      expect(decision.maximum_evidence_distance).to be <= 1e-6 + decision.numeric_tolerance
      verify_monotone_evidence(decision)
    end

    it 'is symmetric when the input order is exchanged' do
      first_result = dense_line.continuous_frechet_distance(shifted_samples_line)
      second_result = shifted_samples_line.continuous_frechet_distance(dense_line)

      expect(second_result.distance).to be_within(1e-12).of(first_result.distance)
      expect(second_result.error_bound).to be_within(1e-12).of(first_result.error_bound)
      expect(shifted_samples_line).to be_continuous_frechet_within_threshold(
        dense_line, 1e-6
      )
    end

    it 'accepts evidence when the threshold exactly touches the free-space boundary' do
      first = GeoRuby::SimpleFeatures::LineString.from_coordinates([[0, 0], [10, 0]])
      second = GeoRuby::SimpleFeatures::LineString.from_coordinates([[0, 1], [10, 1]])

      expect(first.continuous_frechet_within_threshold?(second, 1.0)).to be(true)
      expect(first.continuous_frechet_within_threshold?(second, 0.999)).to be(false)

      decision = first.continuous_frechet_threshold_decision(second, 1.0)
      verify_monotone_evidence(decision)
      expect(decision.maximum_evidence_distance).to be_within(1e-11).of(1.0)
    end

    it 'supports zero-length segments' do
      first = GeoRuby::SimpleFeatures::LineString.from_coordinates([
        [0, 0], [0, 0], [1, 0]
      ])
      second = GeoRuby::SimpleFeatures::LineString.from_coordinates([[0, 0], [1, 0]])

      result = first.continuous_frechet_distance(second)
      expect(result.distance).to be_within(result.error_bound).of(0.0)
      decision = first.continuous_frechet_threshold_decision(second, 0.0)
      verify_monotone_evidence(decision)
    end

    it 'supports a single-point trajectory' do
      point_line = GeoRuby::SimpleFeatures::LineString.from_coordinates([[1, 0]])
      line = GeoRuby::SimpleFeatures::LineString.from_coordinates([
        [0, 0], [1, 0], [3, 0]
      ])

      result = point_line.continuous_frechet_distance(line)
      expect(result.distance).to eq(2.0)
      expect(result.error_bound).to be < 1e-9
      expect(point_line.continuous_frechet_within_threshold?(line, 2.0)).to be(true)
      expect(line.continuous_frechet_within_threshold?(point_line, 2.0)).to be(true)
      expect(point_line.continuous_frechet_within_threshold?(line, 1.999)).to be(false)

      decision = point_line.continuous_frechet_threshold_decision(line, 2.0)
      expect(decision.evidence.map { |step| [step[:u], step[:v]] }).to eq(
        [[0.0, 0.0], [0.0, 0.5], [0.0, 1.0]]
      )
    end

    it 'compares closed trajectories from their supplied first point' do
      first = GeoRuby::SimpleFeatures::LineString.from_coordinates([
        [0, 0], [10, 0], [10, 10], [0, 10], [0, 0]
      ])
      shifted = GeoRuby::SimpleFeatures::LineString.from_coordinates([
        [0, 1], [10, 1], [10, 11], [0, 11], [0, 1]
      ])
      rotated = GeoRuby::SimpleFeatures::LineString.from_coordinates([
        [10, 0], [10, 10], [0, 10], [0, 0], [10, 0]
      ])

      result = first.continuous_frechet_distance(shifted)
      expect(result.distance).to be_within(result.error_bound).of(1.0)
      expect(first.continuous_frechet_within_threshold?(shifted, 1.0)).to be(true)
      expect(first.continuous_frechet_within_threshold?(shifted, 0.999)).to be(false)
      expect(first.continuous_frechet_within_threshold?(rotated, 9.0)).to be(false)
      expect(first.continuous_frechet_within_threshold?(rotated, 10.0)).to be(true)
    end

    it 'rejects planar API calls for three-dimensional line strings' do
      first = GeoRuby::SimpleFeatures::LineString.from_coordinates([[0, 0, 0]], 4326, true)
      second = GeoRuby::SimpleFeatures::LineString.from_coordinates([[0, 0]])

      expect { first.continuous_frechet_distance(second) }
        .to raise_error(ArgumentError, 'continuous Fréchet distance is planar only')
    end
  end

  describe 'Simplify' do

    let(:line) { GeoRuby::SimpleFeatures::LineString.from_coordinates([[6, 0], [4, 1], [3, 4], [4, 6], [5, 8], [5, 9], [4, 10], [6, 15]], 4326) }

    it 'should simplify a simple linestring' do
      line =  GeoRuby::SimpleFeatures::LineString.from_coordinates([[12, 15], [14, 17], [17, 20]], 4326)
      expect(line.simplify.points.size).to eq(2)
    end

    it 'should simplify a harder linestring' do
      expect(line.simplify(6).points.size).to eq(6)
      expect(line.simplify(9).size).to eq(4)
      expect(line.simplify(10).size).to eq(3)
    end

    it 'should barely touch it' do
      expect(line.simplify(1).points.size).to eq(7)
    end

    it 'should simplify to five points' do
      expect(line.simplify(7).points.size).to eq(5)
    end

    it 'should flatten it' do
      expect(line.simplify(11).points.size).to eq(2)
    end

    it 'should be the first and last in a flatten' do
      l = line.simplify(11)
      expect(l[0]).to be_a_point(6, 0)
      expect(l[1]).to be_a_point(6, 15)
    end

  end
end
