require File.expand_path(File.dirname(__FILE__) + '/../../spec_helper')

describe GeoRuby::SimpleFeatures::Frechet do
  let(:frechet) { GeoRuby::SimpleFeatures::Frechet }

  def line(coords)
    GeoRuby::SimpleFeatures::LineString.from_coordinates(coords)
  end

  describe 'continuous vs vertex-based distance' do
    # Tent over a straight segment: the apex vertex (5, 1) sits ~5.1
    # away from every vertex of the other curve, so any vertex-to-vertex
    # (discrete) coupling costs sqrt(26) =~ 5.099, while the continuous
    # Frechet distance is exactly 1 (pair the apex with (5, 0)).
    let(:tent) { line([[0, 0], [5, 1], [10, 0]]) }
    let(:flat) { line([[0, 0], [10, 0]]) }

    it 'is much smaller than the discrete vertex distance' do
      result = tent.frechet_distance(flat)
      expect(result.distance).to be_within(1e-6).of(1.0)
      expect(Math.sqrt(26)).to be > 5 # discrete vertex coupling cost
    end

    it 'stays small when the inputs are swapped' do
      result = flat.frechet_distance(tent)
      expect(result.distance).to be_within(1e-6).of(1.0)
    end

    it 'accepts a threshold exactly touching the free-space boundary' do
      expect(tent.frechet_within?(flat, 1.0)).to be(true)
      expect(flat.frechet_within?(tent, 1.0)).to be(true)
    end

    it 'rejects a threshold just below the distance' do
      expect(tent.frechet_within?(flat, 0.999)).to be(false)
      expect(flat.frechet_within?(tent, 0.999)).to be(false)
    end

    it 'returns a verifiable monotone matching witness' do
      decision = tent.frechet_decision(flat, 1.0)
      expect(decision.within?).to be(true)
      expect(decision.witness.first).to eq([0.0, 0.0])
      expect(decision.witness.last).to eq([2.0, 1.0])
      expect(frechet.valid_witness?(tent.send(:frechet_xy), flat.send(:frechet_xy),
                                    decision.witness, 1.0)).to be(true)
    end
  end

  describe 'distance result' do
    it 'is zero for identical curves' do
      curve = line([[0, 0], [3, 4], [6, 0]])
      result = curve.frechet_distance(line([[0, 0], [3, 4], [6, 0]]))
      expect(result.upper_bound).to be < 1e-6
    end

    it 'equals the constant offset of parallel curves' do
      result = line([[0, 0], [10, 0]]).frechet_distance(line([[0, 3], [10, 3]]))
      expect(result.distance).to be_within(1e-6).of(3.0)
    end

    it 'carries an error bound within the requested tolerance' do
      result = line([[0, 0], [5, 1], [10, 0]]).frechet_distance(line([[0, 0], [10, 0]]), 1e-9)
      expect(result.lower_bound).to be <= 1.0
      expect(result.upper_bound).to be >= 1.0
      expect(result.error_bound).to be <= 1e-9 * [1.0, result.lower_bound].max + 1e-12
    end

    it 'is symmetric up to the error bound' do
      a = line([[0, 0], [2, 3], [5, 1], [8, 4]])
      b = line([[0, 1], [3, 0], [6, 2], [8, 0]])
      ab = a.frechet_distance(b)
      ba = b.frechet_distance(a)
      expect(ab.lower_bound).to be <= ba.upper_bound
      expect(ba.lower_bound).to be <= ab.upper_bound
      expect((ab.distance - ba.distance).abs).to be <= ab.error_bound + ba.error_bound
    end
  end

  describe 'degenerate inputs' do
    it 'handles zero-length segments' do
      with_dup = line([[0, 0], [0, 0], [2, 0], [2, 0]])
      plain = line([[0, 0], [2, 0]])
      expect(with_dup.frechet_distance(plain).upper_bound).to be < 1e-6
      expect(with_dup.frechet_within?(plain, 0.0)).to be(true)
    end

    it 'handles a single-point trajectory against a polyline' do
      point_line = line([[1, 1]])
      other = line([[0, 0], [2, 0]])
      result = point_line.frechet_distance(other)
      expect(result.distance).to be_within(1e-12).of(Math.sqrt(2))
      expect(result.error_bound).to eq(0.0)
    end

    it 'handles a single-point trajectory on both sides' do
      expect(line([[1, 1]]).frechet_distance(line([[4, 5]])).distance).to eq(5.0)
      expect(line([[1, 1]]).frechet_within?(line([[4, 5]]), 5.0)).to be(true)
      expect(line([[1, 1]]).frechet_within?(line([[4, 5]]), 4.999)).to be(false)
    end

    it 'provides a witness for single-point trajectories' do
      decision = line([[1, 1]]).frechet_decision(line([[0, 0], [2, 0]]), Math.sqrt(2))
      expect(decision.within?).to be(true)
      expect(frechet.valid_witness?([[1.0, 1.0]], [[0.0, 0.0], [2.0, 0.0]],
                                    decision.witness, Math.sqrt(2))).to be(true)
    end

    it 'raises on empty line strings' do
      expect { line([[0, 0]]).frechet_distance(GeoRuby::SimpleFeatures::LineString.new) }
        .to raise_error(ArgumentError)
    end
  end

  describe 'closed trajectories' do
    let(:square) { line([[0, 0], [1, 0], [1, 1], [0, 1], [0, 0]]) }
    let(:rotated_square) { line([[1, 0], [1, 1], [0, 1], [0, 0], [1, 0]]) }

    it 'compares them as given, from their own start point' do
      # same geometric loop, but the given start points differ by 1
      result = square.frechet_distance(rotated_square)
      expect(result.lower_bound).to be > 0.99
    end

    it 'is zero for an identical closed trajectory' do
      expect(square.frechet_distance(square).upper_bound).to be < 1e-6
    end
  end

  describe 'witness properties' do
    it 'is monotone and stays within epsilon on a jagged example' do
      a = line([[0, 0], [2, 3], [5, 1], [8, 4], [10, 0]])
      b = line([[0, 1], [3, 0], [6, 2], [8, 0], [10, 1]])
      eps = a.frechet_distance(b).upper_bound
      decision = a.frechet_decision(b, eps)
      expect(decision.within?).to be(true)
      witness = decision.witness
      witness.each_cons(2) do |u, v|
        expect(v[0]).to be >= u[0] - 1e-12
        expect(v[1]).to be >= u[1] - 1e-12
      end
      expect(frechet.valid_witness?(a.send(:frechet_xy), b.send(:frechet_xy),
                                    witness, eps)).to be(true)
    end

    it 'rejects a tampered witness' do
      a = [[0.0, 0.0], [5.0, 1.0], [10.0, 0.0]]
      b = [[0.0, 0.0], [10.0, 0.0]]
      expect(frechet.valid_witness?(a, b, [[0.0, 0.0], [1.0, 0.0], [2.0, 1.0]], 1.0))
        .to be(false)
    end
  end

  describe 'scaling' do
    # decision cost is O(m*n) cell propagations, not an enumeration of
    # vertex alignments; a 120x120 diagram decides quickly
    let(:wave_a) { line((0..120).map { |i| [i * 0.5, Math.sin(i * 0.3)] }) }
    let(:wave_b) { line((0..120).map { |i| [i * 0.5, Math.sin(i * 0.3) + 0.1] }) }

    it 'decides a large diagram' do
      expect(wave_a.frechet_within?(wave_b, 0.2)).to be(true)
      expect(wave_a.frechet_within?(wave_b, 0.05)).to be(false)
    end

    it 'computes a bounded distance' do
      result = wave_a.frechet_distance(wave_b)
      expect(result.distance).to be_within(result.error_bound + 1e-9).of(0.1)
    end
  end
end

