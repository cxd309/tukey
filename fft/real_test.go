package fft

import (
	"math/rand/v2"
	"slices"
	"strconv"
	"testing"

	"github.com/cxd309/tukey/internal/reference"
	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

// randomReal is a seeded random real signal of length n
func randomReal(n int) (x []float64) {
	rng := rand.New(rand.NewPCG(uint64(n), 1))
	x = make([]float64, n)
	for i := range x {
		x[i] = rng.NormFloat64()
	}
	return
}

// the real transforms undo each other,
// for even and odd lengths and through each algorithm:
// radix-2, mixed radix, Bluestein, and Bluestein inside mixed radix
func TestRealRoundTrip(t *testing.T) {
	for _, n := range []int{1, 2, 4, 6, 10, 1024, 1000, 1001, 1009, 2018} {
		t.Run(strconv.Itoa(n), func(t *testing.T) {
			x := randomReal(n)
			y, err := ForwardReal(x)
			require.NoError(t, err)
			require.Len(t, y, n/2+1)
			back, err := InverseReal(y, n)
			require.NoError(t, err)
			reference.AssertClose(t, "x", back, x, reference.Tolerance{Rel: 1e-13, Scaled: 1e-13})
		})
	}
}

// ForwardReal is the first half of Forward on the same signal
func TestForwardRealMatchesForward(t *testing.T) {
	x := randomReal(1000)
	full, err := Forward(complexified(x))
	require.NoError(t, err)
	half, err := ForwardReal(x)
	require.NoError(t, err)
	reference.AssertClose(t, "y", half, full[:501], reference.Tolerance{Rel: 1e-14, Scaled: 1e-14})
}

// DC, and Nyquist for even n, are exactly real, even through Bluestein
func TestForwardRealDCAndNyquistAreReal(t *testing.T) {
	for _, n := range []int{1009, 2018} { // Bluestein, and Bluestein inside mixed radix
		y, err := ForwardReal(randomReal(n))
		require.NoError(t, err)
		assert.Zero(t, imag(y[0]), "DC, n=%d", n)
		if n%2 == 0 {
			assert.Zero(t, imag(y[n/2]), "Nyquist, n=%d", n)
		}
	}
}

// the imaginary parts of DC and Nyquist can't affect a real signal, so they're ignored,
// and bins beyond n/2 + 1 are dropped
func TestInverseRealIgnoresWhatARealSignalCannotHave(t *testing.T) {
	const n = 8
	y, err := ForwardReal(randomReal(n))
	require.NoError(t, err)
	want, err := InverseReal(y, n)
	require.NoError(t, err)

	altered := slices.Clone(y)
	altered[0] += 5i
	altered[n/2] -= 3i
	altered = append(altered, 7+7i, 1i) // extra bins beyond n/2 + 1
	got, err := InverseReal(altered, n)
	require.NoError(t, err)
	assert.Equal(t, want, got)
}

func TestRealTransformsRejectInvalidLengths(t *testing.T) {
	_, err := ForwardReal(nil)
	assert.ErrorIs(t, err, ErrInvalidLength)
	_, err = InverseReal([]complex128{1}, 0)
	assert.ErrorIs(t, err, ErrInvalidLength)
}

// the packed path against the complex one, at even lengths whose half is each
// algorithm: 1 (n=2), odd 3 (6), radix-2, smooth, and Bluestein (2·1009)
func TestForwardRealPackedMatchesComplex(t *testing.T) {
	for _, n := range []int{2, 4, 6, 10, 1024, 1000, 2018} {
		t.Run(strconv.Itoa(n), func(t *testing.T) {
			x := randomReal(n)
			want := transform(complexified(x))[:n/2+1]
			reference.AssertClose(t, "y", forwardRealPacked(x), want, reference.Tolerance{Rel: 1e-14, Scaled: 1e-14})
		})
	}
}
