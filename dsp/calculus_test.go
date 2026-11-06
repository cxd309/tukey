package dsp

import (
	"math"
	"testing"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

// invalidSpacings are sample spacings that aren't positive and finite.
var invalidSpacings = []float64{0, -1, math.NaN(), math.Inf(1), math.Inf(-1)}

func TestDiff(t *testing.T) {
	assert.Equal(t, []float64{3, 5, 7}, Diff([]float64{1, 4, 9, 16}))
	assert.Equal(t, []float64{-2}, Diff([]float64{3, 1}))
}

// Diff of Diff is the second difference: the squares' is constant.
func TestDiffTwice(t *testing.T) {
	assert.Equal(t, []float64{2, 2}, Diff(Diff([]float64{1, 4, 9, 16})))
}

// Fewer than 2 samples have no differences: an empty slice, not nil, so the
// result can always be ranged over, appended to or passed on.
func TestDiffShortSignals(t *testing.T) {
	for _, x := range [][]float64{nil, {}, {5}} {
		d := Diff(x)
		assert.NotNil(t, d)
		assert.Empty(t, d)
	}
}

func TestGrad(t *testing.T) {
	g, err := Grad([]float64{1, 4, 9, 16}, 1)
	require.NoError(t, err)
	assert.Equal(t, []float64{3, 4, 6, 7}, g, "one-sided at the ends, central inside")
}

// A line's gradient is its slope everywhere, ends included. Every sample and
// difference here is exact.
func TestGradOfALineIsItsSlope(t *testing.T) {
	const dx, slope = 0.5, 3.0
	x := make([]float64, 9)
	for i := range x {
		x[i] = 2 + slope*dx*float64(i)
	}
	g, err := Grad(x, dx)
	require.NoError(t, err)
	for i, v := range g {
		assert.Equal(t, slope, v, "g[%d]", i)
	}
}

// With 2 samples both ends are one-sided, the same slope.
func TestGradTwoSamples(t *testing.T) {
	g, err := Grad([]float64{1, 3}, 2)
	require.NoError(t, err)
	assert.Equal(t, []float64{1, 1}, g)
}

func TestGradRejectsInvalidInput(t *testing.T) {
	for _, x := range [][]float64{nil, {}, {5}} {
		_, err := Grad(x, 1)
		assert.ErrorIs(t, err, ErrInvalidLength, "%d samples", len(x))
	}
	for _, dx := range invalidSpacings {
		_, err := Grad([]float64{1, 2, 3}, dx)
		assert.ErrorIs(t, err, ErrInvalidParameter, "dx %g", dx)
	}
}

func TestCumSum(t *testing.T) {
	assert.Equal(t, []float64{1, 3, 6, 10}, CumSum([]float64{1, 2, 3, 4}))
	assert.Equal(t, []float64{-7}, CumSum([]float64{-7}))
}

func TestCumSumEmpty(t *testing.T) {
	for _, x := range [][]float64{nil, {}} {
		s := CumSum(x)
		assert.NotNil(t, s)
		assert.Empty(t, s)
	}
}

// CumSum is compensated as Sum is, so it ends at Sum's total.
func TestCumSumEndsAtSum(t *testing.T) {
	x := mixedMagnitudes()
	s := CumSum(x)
	assert.Equal(t, 20.0, s[len(s)-1])
	assert.Equal(t, 1e8, s[0])
}

func TestTrapz(t *testing.T) {
	cases := []struct {
		name string
		y    []float64
		dx   float64
		want float64
	}{
		{"three samples", []float64{1, 2, 3}, 0.5, 2},
		{"two samples", []float64{2, 4}, 1, 3},
		{"line, exact", []float64{0, 1, 2, 3, 4}, 1, 8},
		{"one sample", []float64{5}, 1, 0},
		{"empty", []float64{}, 1, 0},
		{"nil", nil, 1, 0},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			integral, err := Trapz(c.y, c.dx)
			require.NoError(t, err)
			assert.Equal(t, c.want, integral)
		})
	}
}

func TestCumTrapz(t *testing.T) {
	cases := []struct {
		name string
		y    []float64
		dx   float64
		want []float64
	}{
		{"three samples", []float64{1, 2, 3}, 0.5, []float64{0, 0.75, 2}},
		{"line", []float64{0, 1, 2, 3, 4}, 1, []float64{0, 0.5, 2, 4.5, 8}},
		{"one sample", []float64{5}, 1, []float64{0}},
		{"empty", []float64{}, 1, []float64{}},
		{"nil", nil, 1, []float64{}},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			integral, err := CumTrapz(c.y, c.dx)
			require.NoError(t, err)
			assert.Equal(t, c.want, integral)
		})
	}
}

// CumTrapz ends at Trapz's integral. The two add different terms, the halves
// separately against the middle samples whole, so they can differ in the last
// bit.
func TestCumTrapzEndsAtTrapz(t *testing.T) {
	y, dx := noise(1000, 1), 0.01
	cumulative, err := CumTrapz(y, dx)
	require.NoError(t, err)
	integral, err := Trapz(y, dx)
	require.NoError(t, err)
	assert.InDelta(t, integral, cumulative[len(cumulative)-1], 1e-15*dx*sumAbs(y))
}

// The spacing is checked whatever y's length, so a bad dx is caught even
// before there's anything to integrate.
func TestTrapzRejectsInvalidSpacing(t *testing.T) {
	for _, y := range [][]float64{nil, {1, 2, 3}} {
		for _, dx := range invalidSpacings {
			_, err := Trapz(y, dx)
			assert.ErrorIs(t, err, ErrInvalidParameter, "Trapz, %d samples, dx %g", len(y), dx)
			_, err = CumTrapz(y, dx)
			assert.ErrorIs(t, err, ErrInvalidParameter, "CumTrapz, %d samples, dx %g", len(y), dx)
		}
	}
}

// sumAbs is sum(|x|), the scale a sum of x rounds at.
func sumAbs(x []float64) (s float64) {
	for _, v := range x {
		s += math.Abs(v)
	}
	return
}
