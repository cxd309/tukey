package dsp

import (
	"math"
	"math/big"
	"testing"

	"github.com/cxd309/tukey/window"
	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

// firBands are bands of each type, with where a scaled design has gain 1
var firBands = []struct {
	band    Band
	scaleAt float64
}{
	{Lowpass(0.3), 0},
	{Highpass(0.3), 1},
	{Bandpass(0.2, 0.5), 0.35},
	{Bandstop(0.2, 0.5), 0},
}

// amplitudeAt is the zero-phase amplitude of taps h at f, Σ h·cos(π·m·f)
func amplitudeAt(h []float64, f float64) (a float64) {
	centre := float64(len(h)-1) / 2
	for k, tap := range h {
		a += tap * math.Cos(math.Pi*(float64(k)-centre)*f)
	}
	return
}

// symmetricHamming is a Hamming window made exactly symmetric, by mirroring
func symmetricHamming(t *testing.T, n int) (w []float64) {
	w, err := window.Hamming(n)
	require.NoError(t, err)
	for k := range n / 2 {
		w[n-1-k] = w[k]
	}
	return
}

// A scaled design's response is exactly 1 at DC, Nyquist or its passband's centre.
func TestFIRWinIsScaledToOne(t *testing.T) {
	for _, c := range firBands {
		for _, n := range []int{1, 31, 101} {
			h, err := FIRWin(symmetricHamming(t, n), c.band)
			require.NoError(t, err)
			assert.InDelta(t, 1, amplitudeAt(h, c.scaleAt), 1e-14, "%v, %d taps", c.band, n)
		}
	}
}

// A symmetric window gives exactly symmetric taps, so exactly linear phase:
// sinCosPi is exactly odd and even in its argument.
func TestFIRWinIsSymmetric(t *testing.T) {
	for _, c := range firBands {
		for _, n := range []int{31, 501} {
			h, err := FIRWin(symmetricHamming(t, n), c.band)
			require.NoError(t, err)
			for k := range n / 2 {
				assert.Equal(t, h[k], h[n-1-k], "%v, %d taps, tap %d", c.band, n, k)
			}
		}
	}
}

// An even number of taps has a zero at Nyquist, whatever the band.
func TestFIRWinEvenLengthIsZeroAtNyquist(t *testing.T) {
	for _, band := range []Band{Lowpass(0.3), Bandpass(0.2, 0.5)} {
		h, err := FIRWin(symmetricHamming(t, 32), band)
		require.NoError(t, err)
		assert.InDelta(t, 0, amplitudeAt(h, 1), 1e-15, "%v", band)
	}
}

// With a rectangular window and no scaling, the taps are the ideal response
// truncated: sin(π·cutoff·m)/(π·m), and the cutoff itself at the centre.
func TestFIRWinRectangularIsTruncatedIdeal(t *testing.T) {
	const n, cutoff = 31, 0.3
	ones := make([]float64, n)
	for k := range ones {
		ones[k] = 1
	}
	h, err := FIRWin(ones, Lowpass(cutoff), NoScale())
	require.NoError(t, err)
	for k, tap := range h {
		m := float64(k - n/2)
		want := cutoff
		if m != 0 {
			want = math.Sin(math.Pi*cutoff*m) / (math.Pi * m)
		}
		assert.InDelta(t, want, tap, 1e-16, "tap %d", k)
	}
}

// One tap is the cutoff unscaled and 1 scaled.
func TestFIRWinOneTap(t *testing.T) {
	h, err := FIRWin([]float64{1}, Lowpass(0.3))
	require.NoError(t, err)
	assert.Equal(t, []float64{1}, h)
	h, err = FIRWin([]float64{1}, Lowpass(0.3), NoScale())
	require.NoError(t, err)
	assert.Equal(t, []float64{0.3}, h)
}

func TestFIRWinLeavesWindowUnchanged(t *testing.T) {
	w := symmetricHamming(t, 31)
	original := append([]float64(nil), w...)
	_, err := FIRWin(w, Bandpass(0.2, 0.5))
	require.NoError(t, err)
	assert.Equal(t, original, w)
}

// sinCosPi reduces exactly, so whole and half turns are exact however far
// from 0, and quarter turns are √2/2 to an ULP, as math.Sincos(π/4) is.
func TestSinCosPi(t *testing.T) {
	cases := []struct{ x, sin, cos float64 }{
		{0, 0, 1},
		{0.5, 1, 0},
		{1, 0, -1},
		{-1.5, 1, 0},
		{1000.25, math.Sqrt2 / 2, math.Sqrt2 / 2},
		{-250.5, -1, 0},
		{1e6 + 0.75, math.Sqrt2 / 2, -math.Sqrt2 / 2},
	}
	for _, c := range cases {
		sin, cos := sinCosPi(c.x)
		assert.InDelta(t, c.sin, sin, 2e-16, "sin(π·%g)", c.x)
		assert.InDelta(t, c.cos, cos, 2e-16, "cos(π·%g)", c.x)
	}
}

// sinCosPiProduct adds back a·b's rounding error. The reference reduces the
// exact product mod 2 in big.Float, where it's exact, so only the final
// rounding remains; sin(π·a·b) without the correction drifts by 1e-14 here.
func TestSinCosPiProductKeepsDigits(t *testing.T) {
	const a, b = 0.3, 250.5
	exact := new(big.Float).SetPrec(200).Mul(big.NewFloat(a), big.NewFloat(b))
	turns, _ := new(big.Float).Quo(exact, big.NewFloat(2)).Int(nil)
	reduced := exact.Sub(exact, new(big.Float).SetPrec(200).Mul(new(big.Float).SetInt(turns), big.NewFloat(2)))
	r, _ := reduced.Float64()
	wantSin, wantCos := math.Sincos(math.Pi * r)

	sin, cos := sinCosPiProduct(a, b)
	assert.InDelta(t, wantSin, sin, 2e-16)
	assert.InDelta(t, wantCos, cos, 2e-16)
	naive := math.Sin(math.Pi * a * b)
	assert.Greater(t, math.Abs(naive-wantSin), 1e-15, "the naive form should drift here, or the case tests nothing")
}

func TestFIRWinRejectsInvalidInput(t *testing.T) {
	_, err := FIRWin(nil, Lowpass(0.3))
	assert.ErrorIs(t, err, ErrInvalidLength)
	for _, band := range []Band{Highpass(0.3), Bandstop(0.2, 0.5)} {
		_, err = FIRWin(make([]float64, 32), band)
		assert.ErrorIs(t, err, ErrInvalidLength, "%v with 32 taps", band)
	}
	_, err = FIRWin(make([]float64, 31), Band{})
	assert.ErrorIs(t, err, ErrInvalidBand)
	_, err = FIRWin(make([]float64, 31), Lowpass(1.2))
	assert.ErrorIs(t, err, ErrInvalidBand)
}
