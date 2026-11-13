package dsp

import (
	"math"
	"testing"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

// constantDelay asserts every point of delay is exactly want
func constantDelay(t *testing.T, want float64, delay []float64, what string) {
	t.Helper()
	for i, d := range delay {
		assert.Equal(t, want, d, "%s at point %d", what, i)
	}
}

// A linear-phase FIR's delay is exactly (N − 1)/2 at every frequency, its
// stopband zeros on the unit circle included.
func TestGroupDelayOfLinearPhaseFIR(t *testing.T) {
	for _, taps := range []int{31, 32} {
		h, err := FIRWin(symmetricHamming(t, taps), Lowpass(0.3))
		require.NoError(t, err)
		_, delay, err := GroupDelay(h, []float64{1}, 512)
		require.NoError(t, err)
		constantDelay(t, float64(taps-1)/2, delay, "FIR")
	}
}

// At a zero on the unit circle the phase jumps by π and the delay is its
// limit from either side. A difference [1, −1] has its zero at DC, and a
// 4-tap moving average at π/2, a grid point when n = 512: both are
// antisymmetric or symmetric, so the limit is their constant delay.
func TestGroupDelayLimitAtSymmetricZeros(t *testing.T) {
	_, delay, err := GroupDelay([]float64{1, -1}, []float64{1}, 512)
	require.NoError(t, err)
	constantDelay(t, 0.5, delay, "difference")

	_, delay, err = GroupDelay([]float64{0.25, 0.25, 0.25, 0.25}, []float64{1}, 512)
	require.NoError(t, err)
	constantDelay(t, 1.5, delay, "moving average")
}

// A zero on the unit circle in a polynomial that isn't symmetric takes the
// big.Float path's limit, from F's first nonzero derivative. Each zero on the
// circle adds half a sample, and 1 + x/2 adds Re((x/2)/(1 + x/2)) = 1/3 at
// DC, so (1 − x)(1 + x/2) is 5/6 there and (1 − x)²(1 + x/2), a double zero,
// is 4/3.
func TestGroupDelayLimitAtZerosOnTheCircle(t *testing.T) {
	cases := []struct {
		b    []float64
		want float64
	}{
		{[]float64{1, -0.5, -0.5}, 5.0 / 6},
		{[]float64{1, -1.5, 0, 0.5}, 4.0 / 3},
	}
	for _, c := range cases {
		_, delay, err := GroupDelay(c.b, []float64{1}, 512)
		require.NoError(t, err)
		assert.InDelta(t, c.want, delay[0], 1e-15, "%v at DC", c.b)
		assert.InDelta(t, c.want, delay[1], 1e-4, "%v next to DC, continuous with the limit", c.b)
	}
}

// Each leading zero delays by a whole sample, and trailing zeros change
// nothing.
func TestGroupDelayLeadingAndTrailingZeros(t *testing.T) {
	_, plain, err := GroupDelay([]float64{1, 0.5}, []float64{1, -0.3}, 64)
	require.NoError(t, err)
	_, padded, err := GroupDelay([]float64{0, 0, 1, 0.5, 0}, []float64{1, -0.3, 0}, 64)
	require.NoError(t, err)
	for i := range plain {
		assert.InDelta(t, plain[i]+2, padded[i], 1e-15, "point %d", i)
	}
}

// A denominator (1 − r·x)^8 with r = 63/64 is exactly representable, and its
// delay is 8 times one factor's, −Re(r·x/(1 − r·x)), well conditioned. Near
// DC its own condition number is about 1e17, past delayConditionLimit, so
// those points take the big.Float path; float64 there would be noise. The
// delay crosses 0 near point 29, so it's measured against its largest value.
func TestGroupDelayBigFloatFallback(t *testing.T) {
	const r = 63.0 / 64
	a := []float64{1}
	for range 8 {
		next := make([]float64, len(a)+1)
		for k, c := range a {
			next[k] += c
			next[k+1] -= r * c
		}
		a = next
	}
	w, delay, err := GroupDelay([]float64{1}, a, 512)
	require.NoError(t, err)
	want := make([]float64, len(w))
	var scale float64
	for i, wi := range w {
		x := complex(math.Cos(wi), -math.Sin(wi))
		want[i] = -8 * real(-complex(r, 0)*x/(1-complex(r, 0)*x))
		scale = max(scale, math.Abs(want[i]))
	}
	for i := range w {
		assert.InDelta(t, want[i], delay[i], 1e-14*scale, "point %d", i)
	}
}

// The big.Float path agrees with float64 where float64 is accurate.
func TestBigPolynomialDelayMatchesFloat(t *testing.T) {
	g := delayGrid{64}
	c := []float64{1, -0.4, 0.3, 0.1}
	for i := range g.n {
		assert.InDelta(t, g.polynomialDelay(c, i), g.bigPolynomialDelay(c, i), 1e-15, "point %d", i)
	}
}

// bigSinCosPi agrees with sinCosPi on every quadrant, within sinCosPi's own
// few ULP: π·y and math.Sincos each round, where bigSinCosPi rounds once.
func TestBigSinCosPi(t *testing.T) {
	g := delayGrid{7}
	for num := range 4 * g.n {
		bigSin, bigCos := g.bigSinCosPi(num)
		s, _ := bigSin.Float64()
		c, _ := bigCos.Float64()
		wantSin, wantCos := sinCosPi(float64(num) / float64(2*g.n))
		assert.InDelta(t, wantSin, s, 5*ulp(wantSin), "sin, num %d", num)
		assert.InDelta(t, wantCos, c, 5*ulp(wantCos), "cos, num %d", num)
	}
}

// ulp is the gap between x and the next float64 away from 0
func ulp(x float64) float64 {
	x = math.Abs(x)
	return math.Nextafter(x, math.Inf(1)) - x
}

// SOSGroupDelay of a filter's sections is GroupDelay of its (b, a), for a
// filter whose (b, a) is well conditioned.
func TestSOSGroupDelayMatchesGroupDelay(t *testing.T) {
	f, err := Butter(4, Bandpass(0.2, 0.5))
	require.NoError(t, err)
	b, a, err := f.BA()
	require.NoError(t, err)
	sos, err := f.SOS()
	require.NoError(t, err)
	_, want, err := GroupDelay(b, a, 256)
	require.NoError(t, err)
	_, got, err := SOSGroupDelay(sos, 256)
	require.NoError(t, err)
	assert.InDeltaSlice(t, want, got, 1e-12)
}

// The frequencies are Freqz's, and n = 0 gives empty, non-nil results.
func TestGroupDelayGrid(t *testing.T) {
	b, a := []float64{1, 0.5}, []float64{1, -0.3}
	wantW, _, err := Freqz(b, a, 100)
	require.NoError(t, err)
	w, _, err := GroupDelay(b, a, 100)
	require.NoError(t, err)
	assert.Equal(t, wantW, w)

	w, delay, err := GroupDelay(b, a, 0)
	require.NoError(t, err)
	assert.NotNil(t, w)
	assert.NotNil(t, delay)
	assert.Empty(t, delay)
}

func TestGroupDelayRejectsInvalidInput(t *testing.T) {
	_, _, err := GroupDelay(nil, []float64{1}, 8)
	assert.ErrorIs(t, err, ErrInvalidCoefficients)
	_, _, err = GroupDelay([]float64{1}, []float64{0, 1}, 8)
	assert.ErrorIs(t, err, ErrInvalidCoefficients)
	_, _, err = GroupDelay([]float64{1}, []float64{1}, -1)
	assert.ErrorIs(t, err, ErrInvalidLength)

	_, _, err = SOSGroupDelay(nil, 8)
	assert.ErrorIs(t, err, ErrInvalidSOS)
	_, _, err = SOSGroupDelay(SOS{{1, 0, 0, 2, 0, 0}}, 8)
	assert.ErrorIs(t, err, ErrInvalidSOS)
	_, _, err = SOSGroupDelay(SOS{{1, 0, 0, 1, 0, 0}}, -1)
	assert.ErrorIs(t, err, ErrInvalidLength)
}
