package window

import (
	"math"
	"testing"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

// The general windows reduce to the named ones exactly:
// they share one code path.
func TestSpecialCasesAreExact(t *testing.T) {
	for n := 0; n <= 20; n++ {
		hann, err := Hann(n)
		require.NoError(t, err)
		hamming, err := Hamming(n)
		require.NoError(t, err)

		w, err := GeneralHamming(n, 0.5)
		require.NoError(t, err)
		assert.Equal(t, hann, w, "GeneralHamming(%d, 0.5)", n)

		w, err = GeneralHamming(n, 0.54)
		require.NoError(t, err)
		assert.Equal(t, hamming, w, "GeneralHamming(%d, 0.54)", n)

		w, err = GeneralCosine(n, []float64{0.5, 0.5})
		require.NoError(t, err)
		assert.Equal(t, hann, w, "GeneralCosine(%d, [0.5, 0.5])", n)
	}
}

func TestHannSmall(t *testing.T) {
	w, err := Hann(3)
	require.NoError(t, err)
	assert.Equal(t, []float64{0, 1, 0}, w)

	w, err = Hann(4, Periodic())
	require.NoError(t, err)
	assert.InDeltaSlice(t, []float64{0, 0.5, 1, 0.5}, w, 1e-16)
}

func TestGeneralCosineNeedsCoefficients(t *testing.T) {
	_, err := GeneralCosine(8, nil)
	assert.ErrorIs(t, err, ErrInvalidParameter)
}

func TestCosineSumsRejectNonFinite(t *testing.T) {
	assertInvalidParameter(t, map[string]func() ([]float64, error){
		"GeneralCosine NaN":   func() ([]float64, error) { return GeneralCosine(8, []float64{0.5, math.NaN()}) },
		"GeneralCosine Inf":   func() ([]float64, error) { return GeneralCosine(8, []float64{math.Inf(1)}) },
		"GeneralHamming NaN":  func() ([]float64, error) { return GeneralHamming(8, math.NaN()) },
		"GeneralHamming -Inf": func() ([]float64, error) { return GeneralHamming(8, math.Inf(-1)) },
	})
}
