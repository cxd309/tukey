package dsp

import (
	"slices"
	"testing"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

// A straight line with exactly representable samples detrends to exactly zero.
// Every intermediate is exact: the mean is 4, the centred times are whole
// numbers, and the fitted slope is exactly 0.5.
func TestDetrendLinearRemovesALineExactly(t *testing.T) {
	x := make([]float64, 9)
	for i := range x {
		x[i] = 2 + 0.5*float64(i)
	}
	y, err := Detrend(x, DetrendLinear)
	require.NoError(t, err)
	assert.Equal(t, make([]float64, 9), y)
}

func TestDetrendConstantRemovesTheMean(t *testing.T) {
	y, err := Detrend([]float64{1, 2, 3, 6}, DetrendConstant)
	require.NoError(t, err)
	assert.Equal(t, []float64{-2, -1, 0, 3}, y)
}

// The zero DetrendMode is DetrendLinear, as SciPy's and MATLAB's default.
func TestDetrendZeroModeIsLinear(t *testing.T) {
	var mode DetrendMode
	assert.Equal(t, DetrendLinear, mode)
}

func TestDetrendShortSignals(t *testing.T) {
	for _, mode := range []DetrendMode{DetrendLinear, DetrendConstant} {
		t.Run(mode.String(), func(t *testing.T) {
			y, err := Detrend(nil, mode)
			require.NoError(t, err)
			assert.Empty(t, y)

			y, err = Detrend([]float64{5}, mode)
			require.NoError(t, err)
			assert.Equal(t, []float64{0}, y, "a single sample is its own mean, with no slope")
		})
	}
}

func TestDetrendRejectsUnknownMode(t *testing.T) {
	_, err := Detrend([]float64{1, 2, 3}, DetrendMode(7))
	assert.ErrorIs(t, err, ErrInvalidMode)
	assert.ErrorContains(t, err, "DetrendMode(7)")
}

func TestDetrendDoesNotModifyInput(t *testing.T) {
	x := []float64{1, 4, 2, 8}
	before := slices.Clone(x)
	for _, mode := range []DetrendMode{DetrendLinear, DetrendConstant} {
		_, err := Detrend(x, mode)
		require.NoError(t, err)
		assert.Equal(t, before, x, mode.String())
	}
}
