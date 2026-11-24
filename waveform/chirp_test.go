package waveform

import (
	"math"
	"testing"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

// An empty t gives an empty, non-nil slice.
func TestChirpEmpty(t *testing.T) {
	x, err := Chirp([]float64{}, 1, 1, 80)
	require.NoError(t, err)
	assert.NotNil(t, x)
	assert.Empty(t, x)
}

// Every sweep starts at cos(phi), and Phase(90) turns the cosine into a negated sine.
func TestChirpPhase(t *testing.T) {
	for _, sweep := range []Option{Linear(), Quadratic(), Logarithmic(), Hyperbolic()} {
		x, err := Chirp([]float64{0}, 1, 1, 80, sweep)
		require.NoError(t, err)
		assert.Equal(t, []float64{1}, x)

		x, err = Chirp([]float64{0, 0.25}, 1, 1, 1, sweep, Phase(90))
		require.NoError(t, err)
		assert.InDeltaSlice(t, []float64{0, -1}, x, 1e-16)
	}
}

// At f0 = f1 every sweep is a constant tone,
// including Logarithmic and Hyperbolic, whose formulas divide by 0 there.
func TestChirpConstant(t *testing.T) {
	times := []float64{0, 0.1, 0.25, 1.3, -2}
	want := make([]float64, len(times))
	for i, ti := range times {
		want[i] = math.Cos(2 * math.Pi * 5 * ti)
	}
	for _, sweep := range []Option{Linear(), Quadratic(), Logarithmic(), Hyperbolic()} {
		x, err := Chirp(times, 5, 1, 5, sweep)
		require.NoError(t, err)
		assert.InDeltaSlice(t, want, x, 1e-15)
	}
}

// A NaN or infinite time gives NaN at that sample only.
func TestChirpNonFiniteTime(t *testing.T) {
	for _, sweep := range []Option{Linear(), Quadratic(), Logarithmic(), Hyperbolic()} {
		x, err := Chirp([]float64{0, math.NaN(), math.Inf(1), math.Inf(-1)}, 1, 1, 2, sweep)
		require.NoError(t, err)
		assert.Equal(t, 1.0, x[0])
		for i := 1; i < len(x); i++ {
			assert.True(t, math.IsNaN(x[i]), "sample %d is %g", i, x[i])
		}
	}
}

func TestChirpInvalidParameter(t *testing.T) {
	cases := []struct {
		name       string
		f0, t1, f1 float64
		opts       []Option
	}{
		{"NaN f0", math.NaN(), 1, 2, nil},
		{"infinite t1", 1, math.Inf(1), 2, nil},
		{"NaN f1", 1, 1, math.NaN(), nil},
		{"infinite phase", 1, 1, 2, []Option{Phase(math.Inf(-1))}},
		{"zero t1", 1, 0, 2, nil},
		{"negative t1", 1, -1, 2, nil},
		{"logarithmic zero f0", 0, 1, 2, []Option{Logarithmic()}},
		{"logarithmic opposite signs", -1, 1, 2, []Option{Logarithmic()}},
		{"hyperbolic zero f1", 1, 1, 0, []Option{Hyperbolic()}},
	}
	for _, c := range cases {
		_, err := Chirp([]float64{0}, c.f0, c.t1, c.f1, c.opts...)
		assert.ErrorIs(t, err, ErrInvalidParameter, c.name)
	}
}

// Logarithmic takes negative frequencies of the same sign, same as SciPy.
func TestChirpLogarithmicNegative(t *testing.T) {
	down, err := Chirp([]float64{0, 0.3, 1}, -1, 1, -80, Logarithmic())
	require.NoError(t, err)
	up, err := Chirp([]float64{0, 0.3, 1}, 1, 1, 80, Logarithmic())
	require.NoError(t, err)
	assert.Equal(t, up, down, "cos is even, so negating every frequency changes nothing")
}

func TestChirpInvalidOption(t *testing.T) {
	_, err := Chirp([]float64{0}, 1, 1, 2, VertexAtT1())
	assert.ErrorIs(t, err, ErrInvalidOption, "VertexAtT1 without Quadratic")
	_, err = Chirp([]float64{0}, 1, 1, 2, Logarithmic(), VertexAtT1())
	assert.ErrorIs(t, err, ErrInvalidOption, "VertexAtT1 with Logarithmic")
}

// The last sweep given wins, same as dsp's method options.
func TestChirpLastSweepWins(t *testing.T) {
	times := []float64{0, 0.3, 0.7}
	want, err := Chirp(times, 1, 1, 80, Quadratic())
	require.NoError(t, err)
	got, err := Chirp(times, 1, 1, 80, Logarithmic(), Quadratic())
	require.NoError(t, err)
	assert.Equal(t, want, got)
}
