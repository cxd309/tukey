package dsp

import (
	"math"
	"math/cmplx"
	"testing"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

// notchCases are (w0, q) from narrow to wide: q = 1000 puts the poles hard
// against the unit circle, and the last three have β ≥ sin(π·w0), where the
// poles are real: a centre near Nyquist, β = 1, and a bandwidth past half
// Nyquist
var notchCases = []struct{ w0, q float64 }{
	{0.3, 10},
	{0.05, 35},
	{0.5, 30},
	{0.3, 1000},
	{0.002, 30},
	{0.95, 5},
	{0.1, 1},
	{0.25, 0.5},
	{0.6, 0.8},
}

// notchDesigns are IIRNotch and IIRPeak, with where each passes and stops
var notchDesigns = []struct {
	name   string
	design func(w0, q float64) (ZPK, error)
	passes func(w0 float64) []float64
	stops  func(w0 float64) []float64
}{
	{"IIRNotch", IIRNotch, func(float64) []float64 { return []float64{0, 1} }, func(w0 float64) []float64 { return []float64{w0} }},
	{"IIRPeak", IIRPeak, func(w0 float64) []float64 { return []float64{w0} }, func(float64) []float64 { return []float64{0, 1} }},
}

// halfPowerEdge is where the gain crosses -3 dB between lo and hi, by
// bisection, given it's monotonic there
func halfPowerEdge(f ZPK, lo, hi float64) (edge float64) {
	above := gainDB(f, lo) > halfPowerDB
	for range 100 {
		mid := (lo + hi) / 2
		if (gainDB(f, mid) > halfPowerDB) == above {
			lo = mid
		} else {
			hi = mid
		}
	}
	return (lo + hi) / 2
}

// A notch passes DC and Nyquist at 0 dB and has a null at w0; a peak passes
// w0 at 0 dB and has nulls at DC and Nyquist.
func TestNotchPassesAndStops(t *testing.T) {
	for _, d := range notchDesigns {
		for _, c := range notchCases {
			f, err := d.design(c.w0, c.q)
			require.NoError(t, err)
			for _, w := range d.passes(c.w0) {
				assert.InDelta(t, 0, gainDB(f, w), gainSlackDB, "%s(%g, %g) at %g", d.name, c.w0, c.q, w)
			}
			for _, w := range d.stops(c.w0) {
				assert.Less(t, gainDB(f, w), -200.0, "%s(%g, %g) at %g", d.name, c.w0, c.q, w)
			}
		}
	}
}

// The -3 dB edges either side of w0 are exactly w0/q apart: Orfanidis's
// definition of the bandwidth, which SciPy follows.
func TestNotchBandwidth(t *testing.T) {
	for _, d := range notchDesigns {
		for _, c := range notchCases {
			f, err := d.design(c.w0, c.q)
			require.NoError(t, err)
			low, high := halfPowerEdge(f, 0, c.w0), halfPowerEdge(f, c.w0, 1)
			assert.InEpsilon(t, c.w0/c.q, high-low, 1e-9, "%s(%g, %g)", d.name, c.w0, c.q)
		}
	}
}

// The poles are inside the unit circle and multiply to (1 − β)/(1 + β): a
// conjugate pair while β < sin(π·w0), and two real poles from there.
func TestNotchPoles(t *testing.T) {
	for _, d := range notchDesigns {
		for _, c := range notchCases {
			f, err := d.design(c.w0, c.q)
			require.NoError(t, err)
			require.Len(t, f.Poles, 2)
			for _, p := range f.Poles {
				assert.Less(t, cmplx.Abs(p), 1.0, "%s(%g, %g)", d.name, c.w0, c.q)
			}
			beta := math.Tan(math.Pi * c.w0 / (2 * c.q))
			if beta < math.Sin(math.Pi*c.w0) {
				assert.Equal(t, cmplx.Conj(f.Poles[0]), f.Poles[1], "%s(%g, %g)", d.name, c.w0, c.q)
			} else {
				assert.Zero(t, imag(f.Poles[0]), "%s(%g, %g)", d.name, c.w0, c.q)
				assert.Zero(t, imag(f.Poles[1]), "%s(%g, %g)", d.name, c.w0, c.q)
			}
			assert.InDelta(t, (1-beta)/(1+beta), real(f.Poles[0]*f.Poles[1]), 1e-15, "%s(%g, %g)", d.name, c.w0, c.q)
		}
	}
}

// A narrow notch's poles sit at radius √((1 − β)/(1 + β)), just inside its
// zeros: at q = 1000, 1 − r is about 5e-4, kept to full precision.
func TestNarrowNotchPoleRadius(t *testing.T) {
	const w0, q = 0.3, 1000.0
	f, err := IIRNotch(w0, q)
	require.NoError(t, err)
	beta := math.Tan(math.Pi * w0 / (2 * q))
	want := math.Sqrt((1 - beta) / (1 + beta))
	for _, p := range f.Poles {
		assert.InEpsilon(t, want, cmplx.Abs(p), 1e-15)
		assert.InEpsilon(t, math.Pi*w0, math.Abs(cmplx.Phase(p)), 1e-6, "the poles are near the zeros' angle")
	}
}

// The peak's gain is β/(1 + β); SciPy's 1 − 1/(1 + β) cancels for a narrow
// bandwidth.
func TestIIRPeakGainKeepsDigits(t *testing.T) {
	const w0, q = 0.002, 30.0
	f, err := IIRPeak(w0, q)
	require.NoError(t, err)
	beta := math.Tan(math.Pi * w0 / (2 * q))
	want := beta / (1 + beta)
	assert.InEpsilon(t, want, f.Gain, 1e-15)
	naive := 1 - 1/(1+beta)
	assert.Greater(t, math.Abs(naive-want)/want, 1e-14, "the naive form should drift here, or the case tests nothing")
}

func TestNotchRejectsInvalidInput(t *testing.T) {
	cases := []struct {
		w0, q float64
		what  string
	}{
		{0, 10, "w0 at DC"},
		{1, 10, "w0 at Nyquist"},
		{-0.1, 10, "w0 negative"},
		{1.5, 10, "w0 above Nyquist"},
		{nan, 10, "w0 NaN"},
		{0.3, 0, "q zero"},
		{0.3, -1, "q negative"},
		{0.3, nan, "q NaN"},
		{0.3, math.Inf(1), "q infinite"},
		{0.5, 0.5, "bandwidth at Nyquist"},
		{0.6, 0.5, "bandwidth past Nyquist"},
	}
	for _, d := range notchDesigns {
		for _, c := range cases {
			_, err := d.design(c.w0, c.q)
			assert.ErrorIs(t, err, ErrInvalidParameter, "%s: %s", d.name, c.what)
		}
	}
}
