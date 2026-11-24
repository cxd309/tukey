package waveform

import (
	"math"
	"testing"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

// An empty t gives an empty, non-nil slice.
func TestSweepPolyEmpty(t *testing.T) {
	x, err := SweepPoly([]float64{}, []float64{1})
	require.NoError(t, err)
	assert.NotNil(t, x)
	assert.Empty(t, x)
}

// A linear polynomial is Chirp's linear sweep, and a quadratic with no t term is its Quadratic,
// both exact in double-double, so they agree sample for sample.
func TestSweepPolyIsChirp(t *testing.T) {
	times := make([]float64, 201)
	for i := range times {
		times[i] = float64(i)/50 - 1 // -1 to 3, before 0 and past t1
	}
	cases := []struct {
		name  string
		poly  []float64
		sweep Option
	}{
		{"linear", []float64{79, 1}, Linear()},
		{"quadratic", []float64{79, 0, 1}, Quadratic()},
	}
	for _, c := range cases {
		want, err := Chirp(times, 1, 1, 80, c.sweep, Phase(30))
		require.NoError(t, err)
		got, err := SweepPoly(times, c.poly, Phase(30))
		require.NoError(t, err)
		assert.Equal(t, want, got, c.name)
	}
}

// A constant polynomial is a tone, which Phase(90) turns into a negated sine.
func TestSweepPolyConstant(t *testing.T) {
	x, err := SweepPoly([]float64{0, 0.125, 0.25, 0.5}, []float64{1}, Phase(90))
	require.NoError(t, err)
	assert.InDeltaSlice(t, []float64{0, -math.Sqrt2 / 2, -1, 0}, x, 1e-16)
}

// A NaN or infinite time gives NaN at that sample only.
func TestSweepPolyNonFiniteTime(t *testing.T) {
	x, err := SweepPoly([]float64{0, math.NaN(), math.Inf(1), math.Inf(-1)}, []float64{3, 2, 1})
	require.NoError(t, err)
	assert.Equal(t, 1.0, x[0])
	for i := 1; i < len(x); i++ {
		assert.True(t, math.IsNaN(x[i]), "sample %d is %g", i, x[i])
	}
}

func TestSweepPolyInvalidParameter(t *testing.T) {
	cases := []struct {
		name string
		poly []float64
		opts []Option
	}{
		{"empty poly", []float64{}, nil},
		{"nil poly", nil, nil},
		{"NaN coefficient", []float64{1, math.NaN()}, nil},
		{"infinite coefficient", []float64{math.Inf(1)}, nil},
		{"NaN phase", []float64{1}, []Option{Phase(math.NaN())}},
	}
	for _, c := range cases {
		_, err := SweepPoly([]float64{0}, c.poly, c.opts...)
		assert.ErrorIs(t, err, ErrInvalidParameter, c.name)
	}
}

// SweepPoly's frequency is its polynomial, so Chirp's sweeps don't apply.
func TestSweepPolyInvalidOption(t *testing.T) {
	for _, opt := range []Option{Linear(), Quadratic(), Logarithmic(), Hyperbolic(), VertexAtT1()} {
		_, err := SweepPoly([]float64{0}, []float64{1}, opt)
		assert.ErrorIs(t, err, ErrInvalidOption)
	}
}

// SweepPoly doesn't modify poly.
func TestSweepPolyKeepsPoly(t *testing.T) {
	poly := []float64{0.025, -0.36, 1.25, 2}
	_, err := SweepPoly([]float64{0, 1, 2}, poly)
	require.NoError(t, err)
	assert.Equal(t, []float64{0.025, -0.36, 1.25, 2}, poly)
}
