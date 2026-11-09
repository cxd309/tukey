package dsp

import (
	"math"
	"math/cmplx"
	"testing"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

// gainSlackDB is how far a designed gain may be from its specification.
const gainSlackDB = 1e-8

// gainDB is f's gain in dB at w, normalised to Nyquist.
func gainDB(f ZPK, w float64) (dB float64) {
	z := cmplx.Exp(complex(0, math.Pi*w))
	h := complex(f.Gain, 0)
	for _, zero := range f.Zeros {
		h *= z - zero
	}
	for _, pole := range f.Poles {
		h /= z - pole
	}
	dB = 20 * math.Log10(cmplx.Abs(h))
	return
}

// chebyBands are bands of each type, with the edges their specification applies at.
var chebyBands = []struct {
	name  string
	band  Band
	edges []float64
}{
	{"lowpass", Lowpass(0.3), []float64{0.3}},
	{"highpass", Highpass(0.3), []float64{0.3}},
	{"bandpass", Bandpass(0.2, 0.5), []float64{0.2, 0.5}},
	{"bandstop", Bandstop(0.2, 0.5), []float64{0.2, 0.5}},
}

// A Cheby1 band edge is where the gain first drops below −ripple, not −3 dB.
func TestCheby1GainAtBandEdges(t *testing.T) {
	const ripple = 1.5
	for _, c := range chebyBands {
		for order := 1; order <= 6; order++ {
			f, err := Cheby1(order, ripple, c.band)
			require.NoError(t, err)
			for _, edge := range c.edges {
				assert.InDelta(t, -ripple, gainDB(f, edge), gainSlackDB, "%s, order %d, edge %g", c.name, order, edge)
			}
		}
	}
}

// The passband's peak is 0 dB, so an even order's DC gain is −ripple, and an
// odd order's is 0 dB.
func TestCheby1DCGain(t *testing.T) {
	const ripple = 2.0
	for order := 1; order <= 6; order++ {
		f, err := Cheby1(order, ripple, Lowpass(0.3))
		require.NoError(t, err)
		want := 0.0
		if order%2 == 0 {
			want = -ripple
		}
		assert.InDelta(t, want, gainDB(f, 0), gainSlackDB, "order %d", order)
	}
}

// The passband ripples between 0 dB and −ripple, and never leaves that range.
func TestCheby1PassbandStaysWithinRipple(t *testing.T) {
	const ripple, cutoff = 1.0, 0.3
	f, err := Cheby1(5, ripple, Lowpass(cutoff))
	require.NoError(t, err)
	for i := range 500 {
		w := cutoff * float64(i) / 499
		gain := gainDB(f, w)
		assert.LessOrEqual(t, gain, gainSlackDB, "w %g", w)
		assert.GreaterOrEqual(t, gain, -ripple-gainSlackDB, "w %g", w)
	}
}

// Cheby1 has no finite analog zeros, so a lowpass's are all at Nyquist.
func TestCheby1LowpassZerosAreAtNyquist(t *testing.T) {
	f, err := Cheby1(4, 1, Lowpass(0.3))
	require.NoError(t, err)
	assert.Equal(t, []complex128{-1, -1, -1, -1}, f.Zeros)
}

// A Cheby2 band edge is the start of the stopband, where the gain first
// reaches −attenuation.
func TestCheby2GainAtBandEdges(t *testing.T) {
	const attenuation = 40.0
	for _, c := range chebyBands {
		for order := 1; order <= 6; order++ {
			f, err := Cheby2(order, attenuation, c.band)
			require.NoError(t, err)
			for _, edge := range c.edges {
				assert.InDelta(t, -attenuation, gainDB(f, edge), gainSlackDB, "%s, order %d, edge %g", c.name, order, edge)
			}
		}
	}
}

// The stopband stays at least the attenuation down, and DC is at 0 dB,
// whatever the order's parity.
func TestCheby2Stopband(t *testing.T) {
	const attenuation, cutoff = 30.0, 0.3
	for _, order := range []int{4, 5} {
		f, err := Cheby2(order, attenuation, Lowpass(cutoff))
		require.NoError(t, err)
		assert.InDelta(t, 0, gainDB(f, 0), gainSlackDB, "order %d, DC", order)
		for i := range 500 {
			w := cutoff + (1-cutoff)*float64(i)/500
			assert.LessOrEqual(t, gainDB(f, w), -attenuation+gainSlackDB, "order %d, w %g", order, w)
		}
	}
}

// An odd order's middle analog zero is at infinity, which lands at Nyquist.
// An even order has all its zeros at finite frequencies.
func TestCheby2ZeroAtNyquistOnlyForOddOrders(t *testing.T) {
	odd, err := Cheby2(3, 40, Lowpass(0.3))
	require.NoError(t, err)
	assert.Contains(t, odd.Zeros, complex(-1, 0))
	even, err := Cheby2(4, 40, Lowpass(0.3))
	require.NoError(t, err)
	assert.NotContains(t, even.Zeros, complex(-1, 0))
}

func TestChebyPolesAreStable(t *testing.T) {
	bands := []Band{Lowpass(0.001), Lowpass(0.99), Highpass(0.0008), Bandpass(0.3, 0.31), Bandstop(0.01, 0.99)}
	for _, band := range bands {
		for _, order := range []int{1, 2, 3, 8} {
			one, err := Cheby1(order, 1, band)
			require.NoError(t, err)
			two, err := Cheby2(order, 40, band)
			require.NoError(t, err)
			for _, p := range append(one.Poles, two.Poles...) {
				assert.Less(t, cmplx.Abs(p), 1.0, "band %v, order %d", band, order)
			}
		}
	}
}

// 10^x − 1 cancels almost every digit for a tiny ripple, where Expm1 keeps them.
func TestRippleFactorKeepsDigitsForSmallRipple(t *testing.T) {
	const dB = 1e-12
	x := dB * math.Ln10 / 10
	want := math.Sqrt(x + x*x/2)
	assert.InEpsilon(t, want, rippleFactor(dB), 1e-15)
	naive := math.Sqrt(math.Pow(10, dB/10) - 1)
	assert.Greater(t, math.Abs(naive-want)/want, 1e-6, "the naive form should drift here, or the case tests nothing")
}

func TestChebyRejectsInvalidInput(t *testing.T) {
	for _, dB := range []float64{0, -1, nan, math.Inf(1), math.Inf(-1)} {
		_, err := Cheby1(2, dB, Lowpass(0.3))
		assert.ErrorIs(t, err, ErrInvalidParameter, "Cheby1, %g dB", dB)
		_, err = Cheby2(2, dB, Lowpass(0.3))
		assert.ErrorIs(t, err, ErrInvalidParameter, "Cheby2, %g dB", dB)
	}
	for _, order := range []int{0, -1} {
		_, err := Cheby1(order, 1, Lowpass(0.3))
		assert.ErrorIs(t, err, ErrInvalidOrder, "Cheby1, order %d", order)
		_, err = Cheby2(order, 40, Lowpass(0.3))
		assert.ErrorIs(t, err, ErrInvalidOrder, "Cheby2, order %d", order)
	}
	_, err := Cheby1(2, 1, Band{})
	assert.ErrorIs(t, err, ErrInvalidBand)
	_, err = Cheby2(2, 40, Band{})
	assert.ErrorIs(t, err, ErrInvalidBand)
}
