package dsp

import (
	"math"
	"math/cmplx"
	"testing"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

// A tone around an offset has the offset plus and minus its amplitude as its
// analytic envelopes.
func TestEnvelopeAnalyticTone(t *testing.T) {
	x := make([]float64, 16)
	for i := range x {
		x[i] = 5 + 2*math.Cos(2*math.Pi*3*float64(i)/16)
	}
	upper, lower, err := Envelope(x, EnvelopeAnalytic())
	require.NoError(t, err)
	for i := range x {
		assert.InDelta(t, 7, upper[i], 1e-14, "upper[%d]", i)
		assert.InDelta(t, 3, lower[i], 1e-14, "lower[%d]", i)
	}
}

// The zero EnvelopeMethod is EnvelopeAnalytic.
func TestEnvelopeZeroIsAnalytic(t *testing.T) {
	x := []float64{3, 1, -4, 1, 5, -9, 2, 6}
	wantUpper, wantLower, err := Envelope(x, EnvelopeAnalytic())
	require.NoError(t, err)
	upper, lower, err := Envelope(x, EnvelopeMethod{})
	require.NoError(t, err)
	assert.Equal(t, wantUpper, upper)
	assert.Equal(t, wantLower, lower)
}

func TestEnvelopeRMS(t *testing.T) {
	s := math.Sqrt(4.5)
	cases := []struct {
		name   string
		x      []float64
		window int
		want   []float64 // the upper envelope, and the lower negated, as each x has a mean of 0
	}{
		// a window of 1 is |x - mean|
		{"window 1", []float64{-2, 0, 2}, 1, []float64{2, 0, 2}},
		// an even window has its extra sample after, so k covers k and k+1
		{"even window", []float64{0, 0, 3, -3, 0, 0}, 2, []float64{0, s, 3, s, 0, 0}},
		// the window shrinks at the start rather than counting samples before x
		{"shrinks at the ends", []float64{3, -3, 0, 0, 0, 0}, 3, []float64{3, math.Sqrt(6), math.Sqrt(3), 0, 0, 0}},
		{"longer than x", []float64{1, -1, 1, -1}, 100, []float64{1, 1, 1, 1}},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			upper, lower, err := Envelope(c.x, EnvelopeRMS(c.window))
			require.NoError(t, err)
			assert.Equal(t, c.want, upper)
			for i := range c.want {
				assert.Equal(t, -c.want[i], lower[i], "lower[%d]", i)
			}
		})
	}
}

// The envelopes are centred on x's mean, not on 0.
func TestEnvelopeRMSAroundMean(t *testing.T) {
	upper, lower, err := Envelope([]float64{8, 12, 8, 12}, EnvelopeRMS(100))
	require.NoError(t, err)
	assert.Equal(t, []float64{12, 12, 12, 12}, upper)
	assert.Equal(t, []float64{8, 8, 8, 8}, lower)
}

func TestEnvelopeInvalid(t *testing.T) {
	_, _, err := Envelope(nil, EnvelopeAnalytic())
	assert.ErrorIs(t, err, ErrInvalidLength)
	for _, window := range []int{0, -1} {
		_, _, err = Envelope([]float64{1, 2, 3}, EnvelopeRMS(window))
		assert.ErrorIs(t, err, ErrInvalidParameter)
	}
}

// analyticEnvelope is the magnitude of Hilbert's output.
func TestAnalyticEnvelopeIsMagnitude(t *testing.T) {
	x := []float64{3, 1, -4, 1, 5, -9, 2, 6}
	analytic, err := Hilbert(x)
	require.NoError(t, err)
	env, err := analyticEnvelope(x)
	require.NoError(t, err)
	for i, a := range analytic {
		assert.Equal(t, cmplx.Abs(a), env[i], "env[%d]", i)
	}
}
