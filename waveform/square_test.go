package waveform

import (
	"math"
	"testing"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

// An empty t gives an empty, non-nil slice.
func TestSquareEmpty(t *testing.T) {
	x, err := Square([]float64{})
	require.NoError(t, err)
	assert.NotNil(t, x)
	assert.Empty(t, x)
}

// The default duty is 0.5, so Square is high for the first half of each period, and its edges follow SciPy:
// 1 from the start of a period, -1 from duty·2π.
func TestSquareDefaultDuty(t *testing.T) {
	times := []float64{0, math.Pi / 2, math.Pi, 3 * math.Pi / 2, 2 * math.Pi, -math.Pi / 2}
	want := []float64{1, 1, -1, -1, 1, -1}
	x, err := Square(times)
	require.NoError(t, err)
	assert.Equal(t, want, x)

	withDuty, err := Square(times, Duty(0.5))
	require.NoError(t, err)
	assert.Equal(t, x, withDuty)
}

// Duty 0 is never high, and duty 1 always.
func TestSquareDutyLimits(t *testing.T) {
	times := []float64{0, 1, math.Pi, 6, -2}
	x, err := Square(times, Duty(0))
	require.NoError(t, err)
	assert.Equal(t, []float64{-1, -1, -1, -1, -1}, x)

	x, err = Square(times, Duty(1))
	require.NoError(t, err)
	assert.Equal(t, []float64{1, 1, 1, 1, 1}, x)
}

// A NaN or infinite time gives NaN at that sample only, where SciPy and MATLAB give -1.
func TestSquareNonFiniteTime(t *testing.T) {
	x, err := Square([]float64{0, math.NaN(), math.Inf(1), math.Inf(-1)})
	require.NoError(t, err)
	assert.Equal(t, 1.0, x[0])
	for i := 1; i < len(x); i++ {
		assert.True(t, math.IsNaN(x[i]), "sample %d is %g", i, x[i])
	}
}

func TestSquareInvalidDuty(t *testing.T) {
	for _, duty := range []float64{-0.1, 1.1, math.NaN(), math.Inf(1), math.Inf(-1)} {
		_, err := Square([]float64{0}, Duty(duty))
		assert.ErrorIs(t, err, ErrInvalidParameter, "duty %g", duty)
	}
}

func TestSquareInvalidOption(t *testing.T) {
	for _, opt := range []Option{Phase(90), Linear(), VertexAtT1()} {
		_, err := Square([]float64{0}, opt)
		assert.ErrorIs(t, err, ErrInvalidOption)
	}
}

// Duty doesn't apply to Chirp or SweepPoly.
func TestDutyInvalidElsewhere(t *testing.T) {
	_, err := Chirp([]float64{0}, 1, 1, 2, Duty(0.5))
	assert.ErrorIs(t, err, ErrInvalidOption, "Chirp")
	_, err = SweepPoly([]float64{0}, []float64{1}, Duty(0.5))
	assert.ErrorIs(t, err, ErrInvalidOption, "SweepPoly")
}
