package fft

import (
	"fmt"
	"math"
	"math/rand/v2"
	"slices"
	"testing"

	"github.com/cxd309/tukey/internal/reference"
	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

// transformTolerance was measured against scipy.fft:
//   - radix-2 and mixed radix both within 2 ULP, worst 4.1e-16, on every length they handle
//   - Bluestein within 4 ULP, not growing with n up to 4096
//
// It's about 10 times the worst.
var transformTolerance = reference.Tolerance{Rel: 1e-14, Scaled: 1e-14}

// randomSignal is a seeded random complex signal of length n.
func randomSignal(n int) (x []complex128) {
	rng := rand.New(rand.NewPCG(uint64(n), 0))
	x = make([]complex128, n)
	for i := range x {
		x[i] = complex(rng.NormFloat64(), rng.NormFloat64())
	}
	return
}

// checkAgainstSciPy runs transform on every vector in category.
func checkAgainstSciPy(t *testing.T, category string, transform func([]complex128, ...Option) ([]complex128, error)) {
	reference.Run(t, reference.SourceSciPy, category, func(t *testing.T, v transformVector) {
		x, err := v.Params.X.Values()
		require.NoError(t, err)
		want, err := v.Output.Y.Values()
		require.NoError(t, err)
		got, err := transform(x)
		require.NoError(t, err)
		reference.AssertClose(t, "y", got, want, transformTolerance)
	})
}

func TestForwardMatchesSciPy(t *testing.T) { checkAgainstSciPy(t, "fft", Forward) }
func TestInverseMatchesSciPy(t *testing.T) { checkAgainstSciPy(t, "ifft", Inverse) }

// Parseval's theorem:
// energy is the same in both domains, Σ|x|² = (1/n)·Σ|X|².
func TestForwardPreservesEnergy(t *testing.T) {
	x := randomSignal(512)
	y, err := Forward(x)
	require.NoError(t, err)
	timeEnergy, freqEnergy := 0.0, 0.0
	for i := range x {
		timeEnergy += real(x[i])*real(x[i]) + imag(x[i])*imag(x[i])
		freqEnergy += real(y[i])*real(y[i]) + imag(y[i])*imag(y[i])
	}
	assert.InEpsilon(t, timeEnergy, freqEnergy/512, 1e-13)
}

// An impulse at 0 contains every frequency equally, so its transform is all ones, exactly:
// each butterfly only ever adds a value to a product with exactly zero.
func TestForwardImpulseIsAllOnes(t *testing.T) {
	x := make([]complex128, 16)
	x[0] = 1
	y, err := Forward(x)
	require.NoError(t, err)
	assert.Equal(t, slices.Repeat([]complex128{1}, 16), y)
}

func TestForwardRejectsEmpty(t *testing.T) {
	_, err := Forward(nil)
	assert.ErrorIs(t, err, ErrSignalTooShort)
}

func TestForwardDoesNotModifyInput(t *testing.T) {
	x := randomSignal(64)
	before := slices.Clone(x)
	_, err := Forward(x)
	require.NoError(t, err)
	assert.Equal(t, before, x)
}

// dft is the discrete Fourier transform straight from its definition, in O(n²):
//
//	y[k] = Σ x[j]·e^(-2πi·jk/n)
//
// It's far too slow for real use, and is the reference the fast algorithms are checked against.
func dft(x []complex128) (y []complex128) {
	n := len(x)
	y = make([]complex128, n)
	for k := range y {
		var sum complex128
		for j, xj := range x {
			sum += xj * twiddle(j*k%n, n)
		}
		y[k] = sum
	}
	return
}

type transformVector struct {
	reference.Meta
	Params struct {
		X reference.Complex `json:"x"`
	} `json:"params"`
	Output struct {
		Y reference.Complex `json:"y"`
	} `json:"output"`
}

// The oracle must itself agree with SciPy before anything is tested against it.
// Measured: 0 to 2 ULP at small n, growing with n as its O(n) summations accumulate rounding,
// to 18 ULP, 4.0e-15, at n = 4096.
// Tolerance: about 10 times the worst.
func TestDFTMatchesSciPy(t *testing.T) {
	reference.Run(t, reference.SourceSciPy, "fft", func(t *testing.T, v transformVector) {
		x, err := v.Params.X.Values()
		require.NoError(t, err)
		want, err := v.Output.Y.Values()
		require.NoError(t, err)
		reference.AssertClose(t, "y", dft(x), want, reference.Tolerance{Rel: 1e-13, Scaled: 1e-13})
	})
}

// benchmarkLengths covers each algorithm path, at a small and large size.
var benchmarkLengths = []struct {
	path string
	n    int
}{
	{"radix2", 1024},
	{"radix2", 65536},
	{"smooth", 1000},
	{"smooth", 48000}, // 2^7x3x5^3: one second at 48 kHz
	{"prime", 1009},   // a single Bluestein butterfly
	{"prime", 65537},
	{"largePrimeFactor", 4036}, // 4x1009: Bluestein butterflies inside mixed radix
}

func BenchmarkForward(b *testing.B) {
	for _, c := range benchmarkLengths {
		b.Run(fmt.Sprintf("%s/n=%d", c.path, c.n), func(b *testing.B) {
			x := randomSignal(c.n)
			b.ReportAllocs()
			b.ResetTimer()
			for range b.N {
				if _, err := Forward(x); err != nil {
					b.Fatal(err)
				}
			}
			reportPerNLogN(b, c.n)
		})
	}
}

// reportPerNLogN adds the time per operation divided by n·log2(n).
// For an efficient FFT that's about constant across lengths, so it shows up an algorithm that's slow for its size.
func reportPerNLogN(b *testing.B, n int) {
	perOp := float64(b.Elapsed().Nanoseconds()) / float64(b.N)
	b.ReportMetric(perOp/(float64(n)*math.Log2(float64(n))), "ns/nlogn")
}
