package window

import (
	"math"
	"testing"

	"github.com/cxd309/tukey/internal/special"
	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

// beta 0 is Boxcar exactly: every sample is I0(0) / I0(0)
func TestKaiserAtBeta0IsBoxcar(t *testing.T) {
	for n := 0; n <= 20; n++ {
		boxcar, err := Boxcar(n)
		require.NoError(t, err)
		w, err := Kaiser(n, 0)
		require.NoError(t, err)
		assert.Equal(t, boxcar, w, "n = %d", n)
	}
}

// I0 is even, so a negative beta is the same window
func TestKaiserNegativeBetaIsSameWindow(t *testing.T) {
	positive, err := Kaiser(16, 5)
	require.NoError(t, err)
	negative, err := Kaiser(16, -5)
	require.NoError(t, err)
	assert.Equal(t, positive, negative)
}

// the ends of a symmetric Kaiser window are I0(0) / I0(beta) = 1 / I0(beta)
func TestKaiserEnds(t *testing.T) {
	w, err := Kaiser(9, 8.6)
	require.NoError(t, err)
	assert.Equal(t, 1/special.I0(8.6), w[0])
	assert.Equal(t, 1.0, w[4], "the centre is I0(beta) / I0(beta)")
}

// the Princen–Bradley condition, which is what KBD is for:
// w[k]² + w[k+n/2]² = 1, so overlapped MDCT frames reconstruct perfectly
func TestKaiserBesselDerivedPrincenBradley(t *testing.T) {
	for _, n := range []int{2, 8, 64, 256} {
		w, err := KaiserBesselDerived(n, 4)
		require.NoError(t, err)
		for k := range n / 2 {
			assert.InDelta(t, 1, w[k]*w[k]+w[k+n/2]*w[k+n/2], 1e-15, "n = %d, k = %d", n, k)
		}
	}
}

func TestKaiserBesselDerivedLengths(t *testing.T) {
	w, err := KaiserBesselDerived(0, 4)
	require.NoError(t, err)
	assert.Empty(t, w)

	for _, n := range []int{-2, 1, 3, 9} {
		_, err := KaiserBesselDerived(n, 4)
		assert.ErrorIs(t, err, ErrInvalidLength, "n = %d", n)
	}
}

func TestKaiserRejectsInvalidParameters(t *testing.T) {
	assertInvalidParameter(t, map[string]func() ([]float64, error){
		"Kaiser beta NaN":              func() ([]float64, error) { return Kaiser(8, math.NaN()) },
		"KaiserBesselDerived beta Inf": func() ([]float64, error) { return KaiserBesselDerived(8, math.Inf(1)) },
	})
}
