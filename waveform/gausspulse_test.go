package waveform

import (
	"math"
	"slices"
	"testing"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

// An empty t gives empty, non-nil slices.
func TestGaussPulseEmpty(t *testing.T) {
	pulse, err := GaussPulse([]float64{}, 1000)
	require.NoError(t, err)
	for _, x := range [][]float64{pulse.InPhase, pulse.Quadrature, pulse.Envelope} {
		assert.NotNil(t, x)
		assert.Empty(t, x)
	}
}

// At t = 0 the envelope peaks at 1, with the carrier's cosine at 1 and its sine at 0.
func TestGaussPulsePeak(t *testing.T) {
	pulse, err := GaussPulse([]float64{0}, 1000)
	require.NoError(t, err)
	assert.Equal(t, GaussianPulse{InPhase: []float64{1}, Quadrature: []float64{0}, Envelope: []float64{1}}, pulse)
}

// The defaults are Bandwidth(0.5) and BandwidthLevel(-6), and CutoffLevel(-60) for the cutoff.
func TestGaussPulseDefaults(t *testing.T) {
	times := []float64{-1e-3, -2.5e-4, 0, 3e-4, 1.1e-3}
	want, err := GaussPulse(times, 1000, Bandwidth(0.5), BandwidthLevel(-6))
	require.NoError(t, err)
	got, err := GaussPulse(times, 1000)
	require.NoError(t, err)
	assert.Equal(t, want, got)

	wantTc, err := GaussPulseCutoff(1000, Bandwidth(0.5), BandwidthLevel(-6), CutoffLevel(-60))
	require.NoError(t, err)
	gotTc, err := GaussPulseCutoff(1000)
	require.NoError(t, err)
	assert.Equal(t, wantTc, gotTc)
}

// The envelope at the cutoff is the cutoff level, 10^(level/20), for each shape.
func TestGaussPulseCutoffLevel(t *testing.T) {
	cases := []struct {
		name  string
		shape []Option
		db    float64
	}{
		{"defaults", nil, -60},
		{"ultrasound", []Option{Bandwidth(0.6)}, -40},
		{"narrow", []Option{Bandwidth(0.1), BandwidthLevel(-3)}, -100},
		{"at the bandwidth level", nil, -6},
	}
	for _, c := range cases {
		tc, err := GaussPulseCutoff(50e3, slices.Concat(c.shape, []Option{CutoffLevel(c.db)})...)
		require.NoError(t, err, c.name)
		pulse, err := GaussPulse([]float64{-tc, tc}, 50e3, c.shape...)
		require.NoError(t, err, c.name)
		want := math.Pow(10, c.db/20)
		assert.InEpsilonSlice(t, []float64{want, want}, pulse.Envelope, 1e-14, c.name)
	}
}

// The envelope and in-phase part are even in t, and the quadrature part odd.
func TestGaussPulseSymmetric(t *testing.T) {
	times := []float64{1e-5, 1.3e-4, 2.5e-4, 7e-4, 1.9e-3}
	negated := make([]float64, len(times))
	for i, ti := range times {
		negated[i] = -ti
	}
	pulse, err := GaussPulse(times, 1000)
	require.NoError(t, err)
	mirror, err := GaussPulse(negated, 1000)
	require.NoError(t, err)
	assert.Equal(t, pulse.Envelope, mirror.Envelope)
	assert.InDeltaSlice(t, pulse.InPhase, mirror.InPhase, 1e-16)
	for i := range times {
		assert.InDelta(t, -pulse.Quadrature[i], mirror.Quadrature[i], 1e-16, "sample %d", i)
	}
}

// A NaN time gives NaN at that sample only, and an infinite time 0, the limit of the pulse.
func TestGaussPulseNonFiniteTime(t *testing.T) {
	pulse, err := GaussPulse([]float64{0, math.NaN(), math.Inf(1), math.Inf(-1)}, 1000)
	require.NoError(t, err)
	for _, x := range [][]float64{pulse.InPhase, pulse.Envelope} {
		assert.Equal(t, 1.0, x[0])
		assert.True(t, math.IsNaN(x[1]), "sample 1 is %g", x[1])
		assert.Equal(t, []float64{0, 0}, x[2:])
	}
	assert.True(t, math.IsNaN(pulse.Quadrature[1]), "sample 1 is %g", pulse.Quadrature[1])
	assert.Equal(t, []float64{0, 0}, pulse.Quadrature[2:])
}

// Both functions check fc and the bandwidth the same way.
func TestGaussPulseInvalidParameter(t *testing.T) {
	cases := []struct {
		name string
		fc   float64
		opts []Option
	}{
		{"zero fc", 0, nil},
		{"negative fc", -1000, nil},
		{"NaN fc", math.NaN(), nil},
		{"infinite fc", math.Inf(1), nil},
		{"zero Bandwidth", 1000, []Option{Bandwidth(0)}},
		{"negative Bandwidth", 1000, []Option{Bandwidth(-0.5)}},
		{"NaN Bandwidth", 1000, []Option{Bandwidth(math.NaN())}},
		{"infinite Bandwidth", 1000, []Option{Bandwidth(math.Inf(1))}},
		{"zero BandwidthLevel", 1000, []Option{BandwidthLevel(0)}},
		{"positive BandwidthLevel", 1000, []Option{BandwidthLevel(3)}},
		{"NaN BandwidthLevel", 1000, []Option{BandwidthLevel(math.NaN())}},
		{"infinite BandwidthLevel", 1000, []Option{BandwidthLevel(math.Inf(-1))}},
	}
	for _, c := range cases {
		_, err := GaussPulse([]float64{0}, c.fc, c.opts...)
		assert.ErrorIs(t, err, ErrInvalidParameter, "GaussPulse %s", c.name)
		_, err = GaussPulseCutoff(c.fc, c.opts...)
		assert.ErrorIs(t, err, ErrInvalidParameter, "GaussPulseCutoff %s", c.name)
	}
}

func TestGaussPulseCutoffInvalidLevel(t *testing.T) {
	for _, db := range []float64{0, 6, math.NaN(), math.Inf(-1)} {
		_, err := GaussPulseCutoff(1000, CutoffLevel(db))
		assert.ErrorIs(t, err, ErrInvalidParameter, "CutoffLevel %g", db)
	}
}

// CutoffLevel applies only to GaussPulseCutoff, and neither takes the other waveforms' options.
func TestGaussPulseInvalidOption(t *testing.T) {
	for _, opt := range []Option{CutoffLevel(-60), Duty(0.5), Phase(90), Linear()} {
		_, err := GaussPulse([]float64{0}, 1000, opt)
		assert.ErrorIs(t, err, ErrInvalidOption)
	}
	for _, opt := range []Option{Duty(0.5), Phase(90), Linear()} {
		_, err := GaussPulseCutoff(1000, opt)
		assert.ErrorIs(t, err, ErrInvalidOption)
	}
}

// Bandwidth, BandwidthLevel and CutoffLevel don't apply to the other waveforms.
func TestGaussOptionsInvalidElsewhere(t *testing.T) {
	for _, opt := range []Option{Bandwidth(0.5), BandwidthLevel(-6), CutoffLevel(-60)} {
		_, err := Chirp([]float64{0}, 1, 1, 2, opt)
		assert.ErrorIs(t, err, ErrInvalidOption, "Chirp")
		_, err = SweepPoly([]float64{0}, []float64{1}, opt)
		assert.ErrorIs(t, err, ErrInvalidOption, "SweepPoly")
		_, err = Square([]float64{0}, opt)
		assert.ErrorIs(t, err, ErrInvalidOption, "Square")
		_, err = Sawtooth([]float64{0}, opt)
		assert.ErrorIs(t, err, ErrInvalidOption, "Sawtooth")
	}
}
