package waveform

import (
	"math"
	"testing"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

// An empty t gives an empty, non-nil slice.
func TestSawtoothEmpty(t *testing.T) {
	x, err := Sawtooth([]float64{})
	require.NoError(t, err)
	assert.NotNil(t, x)
	assert.Empty(t, x)
}

// The default duty is 1, a ramp rising from -1 at the start of each period, same as Duty(1).
func TestSawtoothDefaultDuty(t *testing.T) {
	times := []float64{0, math.Pi / 2, math.Pi, -math.Pi / 2, 2 * math.Pi}
	want := []float64{-1, -0.5, 0, 0.5, -1}
	x, err := Sawtooth(times)
	require.NoError(t, err)
	assert.Equal(t, want, x)

	withDuty, err := Sawtooth(times, Duty(1))
	require.NoError(t, err)
	assert.Equal(t, x, withDuty)
}

// Duty(0.5) is a triangle, peaking at π, and Duty(0) a falling ramp.
func TestSawtoothShapes(t *testing.T) {
	x, err := Sawtooth([]float64{0, math.Pi / 2, math.Pi, -math.Pi}, Duty(0.5))
	require.NoError(t, err)
	assert.Equal(t, []float64{-1, 0, 1, 1}, x)

	x, err = Sawtooth([]float64{0, math.Pi / 2, math.Pi}, Duty(0))
	require.NoError(t, err)
	assert.Equal(t, []float64{1, 0.5, 0}, x)
}

// A tiny negative t is just before the end of a period, so a rising ramp is just below 1,
// where SciPy's t % 2π rounds up to 2π and gives NaN.
func TestSawtoothTinyNegative(t *testing.T) {
	x, err := Sawtooth([]float64{-1e-300, -1e-20, -0.0})
	require.NoError(t, err)
	assert.Equal(t, []float64{1, 1, -1}, x)
}

// A NaN or infinite time gives NaN at that sample only.
func TestSawtoothNonFiniteTime(t *testing.T) {
	x, err := Sawtooth([]float64{0, math.NaN(), math.Inf(1), math.Inf(-1)})
	require.NoError(t, err)
	assert.Equal(t, -1.0, x[0])
	for i := 1; i < len(x); i++ {
		assert.True(t, math.IsNaN(x[i]), "sample %d is %g", i, x[i])
	}
}

func TestSawtoothInvalidDuty(t *testing.T) {
	for _, duty := range []float64{-0.1, 1.1, math.NaN(), math.Inf(1), math.Inf(-1)} {
		_, err := Sawtooth([]float64{0}, Duty(duty))
		assert.ErrorIs(t, err, ErrInvalidParameter, "duty %g", duty)
	}
}

func TestSawtoothInvalidOption(t *testing.T) {
	for _, opt := range []Option{Phase(90), Linear(), VertexAtT1()} {
		_, err := Sawtooth([]float64{0}, opt)
		assert.ErrorIs(t, err, ErrInvalidOption)
	}
}
