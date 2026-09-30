package fft

import (
	"math/rand/v2"
	"slices"
	"testing"

	"github.com/cxd309/tukey/internal/reference"
	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

// radix2Tolerance was measured against scipy.fft on every power-of-two vector:
// worst 4.1e-16 (2 ULP) forward and 3.6e-16 inverse, not growing with n up to
// 4096; set ~10x above the worst
var radix2Tolerance = reference.Tolerance{Rel: 1e-14, Scaled: 1e-14}

// randomSignal is a seeded random complex signal of length n
func randomSignal(n int) (x []complex128) {
	rng := rand.New(rand.NewPCG(uint64(n), 0))
	x = make([]complex128, n)
	for i := range x {
		x[i] = complex(rng.NormFloat64(), rng.NormFloat64())
	}
	return
}

// checkAgainstSciPy runs transform on every vector in category, skipping lengths
// the package doesn't support yet
func checkAgainstSciPy(t *testing.T, category string, transform func([]complex128) ([]complex128, error)) {
	reference.Run(t, category, func(t *testing.T, v transformVector) {
		x, err := v.Params.X.Values()
		require.NoError(t, err)
		if !isPowerOfTwo(len(x)) {
			t.Skip("only powers of two until mixed radix and Bluestein")
		}
		want, err := v.Output.Y.Values()
		require.NoError(t, err)
		got, err := transform(x)
		require.NoError(t, err)
		reference.AssertClose(t, "y", got, want, radix2Tolerance)
	})
}

func TestForwardMatchesSciPy(t *testing.T) { checkAgainstSciPy(t, "fft", Forward) }
func TestInverseMatchesSciPy(t *testing.T) { checkAgainstSciPy(t, "ifft", Inverse) }

// the fast transform must agree with the definition, at every power of two
func TestForwardMatchesDFT(t *testing.T) {
	for n := 1; n <= 256; n *= 2 {
		x := randomSignal(n)
		y, err := Forward(x)
		require.NoError(t, err)
		reference.AssertClose(t, "y", y, dft(x), reference.Tolerance{Rel: 1e-13, Scaled: 1e-13})
	}
}

// the inverse undoes the forward transform
func TestInverseRoundTrip(t *testing.T) {
	x := randomSignal(1024)
	y, err := Forward(x)
	require.NoError(t, err)
	back, err := Inverse(y)
	require.NoError(t, err)
	reference.AssertClose(t, "x", back, x, reference.Tolerance{Rel: 1e-13, Scaled: 1e-13})
}

// Parseval's theorem: energy is the same in both domains, Σ|x|² = (1/n)·Σ|X|²
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

// an impulse at 0 contains every frequency equally, so its transform is all ones,
// exactly: each butterfly only ever adds a value to a product with exactly zero
func TestForwardImpulseIsAllOnes(t *testing.T) {
	x := make([]complex128, 16)
	x[0] = 1
	y, err := Forward(x)
	require.NoError(t, err)
	assert.Equal(t, slices.Repeat([]complex128{1}, 16), y)
}

func TestForwardRejectsEmpty(t *testing.T) {
	_, err := Forward(nil)
	assert.ErrorIs(t, err, ErrInvalidLength)
}

func TestForwardDoesNotModifyInput(t *testing.T) {
	x := randomSignal(64)
	before := slices.Clone(x)
	_, err := Forward(x)
	require.NoError(t, err)
	assert.Equal(t, before, x)
}
