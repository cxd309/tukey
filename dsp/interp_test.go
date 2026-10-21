package dsp

import (
	"math"
	"slices"
	"testing"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

// On a sample point Interp1 returns y itself, exactly, including at both ends.
func TestInterp1OnSamplePointsIsExact(t *testing.T) {
	x := []float64{0, 0.1, 0.35, 1, 2.7}
	y := []float64{1.1, -3.3, 0.7, 9.9, 0.01}
	yi, err := Interp1(x, y, x)
	require.NoError(t, err)
	assert.Equal(t, y, yi)
}

func TestInterp1Between(t *testing.T) {
	yi, err := Interp1([]float64{0, 2, 4}, []float64{0, 4, 2}, []float64{1, 3, 0.5})
	require.NoError(t, err)
	assert.Equal(t, []float64{2, 3, 1}, yi)
}

// Outside [x[0], x[len(x)-1]], and for NaN, the result is NaN, as MATLAB's
// interp1 and SciPy's interp1d with bounds_error=False.
func TestInterp1OutsideIsNaN(t *testing.T) {
	yi, err := Interp1([]float64{0, 1, 2}, []float64{5, 6, 7}, []float64{-0.001, 2.001, math.NaN(), math.Inf(1)})
	require.NoError(t, err)
	for i, v := range yi {
		assert.True(t, math.IsNaN(v), "yi[%d] = %g", i, v)
	}
}

func TestInterp1NoQueryPoints(t *testing.T) {
	yi, err := Interp1([]float64{0, 1}, []float64{0, 1}, nil)
	require.NoError(t, err)
	assert.Empty(t, yi)
}

func TestInterp1RejectsInvalidInput(t *testing.T) {
	cases := []struct {
		name    string
		x, y    []float64
		wantErr error
	}{
		{"different lengths", []float64{0, 1, 2}, []float64{0, 1}, ErrInvalidLength},
		{"one sample", []float64{0}, []float64{0}, ErrInvalidLength},
		{"repeated x", []float64{0, 1, 1}, []float64{0, 1, 2}, ErrInvalidSamplePoints},
		{"decreasing x", []float64{2, 1, 0}, []float64{0, 1, 2}, ErrInvalidSamplePoints},
		{"NaN in x", []float64{0, math.NaN(), 2}, []float64{0, 1, 2}, ErrInvalidSamplePoints},
		{"NaN first", []float64{math.NaN(), 1, 2}, []float64{0, 1, 2}, ErrInvalidSamplePoints},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			_, err := Interp1(c.x, c.y, []float64{0.5})
			assert.ErrorIs(t, err, c.wantErr)
		})
	}
}

func TestInterp1DoesNotModifyInputs(t *testing.T) {
	x, y, xi := []float64{0, 1, 2}, []float64{3, 1, 4}, []float64{1.5, 0.5}
	xBefore, yBefore, xiBefore := slices.Clone(x), slices.Clone(y), slices.Clone(xi)

	_, err := Interp1(x, y, xi)
	require.NoError(t, err)

	assert.Equal(t, xBefore, x)
	assert.Equal(t, yBefore, y)
	assert.Equal(t, xiBefore, xi)
}
