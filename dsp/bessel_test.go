package dsp

import (
	"math"
	"math/cmplx"
	"testing"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

var besselNorms = []BesselNorm{BesselNormMagnitude, BesselNormPhase, BesselNormDelay}

// halfPowerDB is the gain at half power, 10·log10(1/2)
var halfPowerDB = 10 * math.Log10(0.5)

// groupDelayAtDC is f's group delay at DC in samples: each pole p adds
// Re(1/(1 − p)) and each zero z takes away Re(1/(1 − z)), the derivative of
// their phase at ω = 0.
func groupDelayAtDC(f ZPK) (delay float64) {
	for _, p := range f.Poles {
		delay += real(1 / (1 - p))
	}
	for _, z := range f.Zeros {
		delay -= real(1 / (1 - z))
	}
	return
}

// By default a band edge is the -3 dB point, as for Butter.
func TestBesselGainAtBandEdges(t *testing.T) {
	for _, c := range chebyBands {
		for order := 1; order <= 8; order++ {
			f, err := Bessel(order, c.band)
			require.NoError(t, err)
			for _, edge := range c.edges {
				assert.InDelta(t, halfPowerDB, gainDB(f, edge), gainSlackDB, "%s, order %d, edge %g", c.name, order, edge)
			}
		}
	}
}

// In every normalisation, a lowpass passes DC and a highpass Nyquist at 0 dB.
func TestBesselPassbandGain(t *testing.T) {
	for _, norm := range besselNorms {
		for order := 1; order <= 8; order++ {
			low, err := Bessel(order, Lowpass(0.3), BesselNormalised(norm))
			require.NoError(t, err)
			assert.InDelta(t, 0, gainDB(low, 0), gainSlackDB, "norm %d, order %d", norm, order)
			high, err := Bessel(order, Highpass(0.3), BesselNormalised(norm))
			require.NoError(t, err)
			assert.InDelta(t, 0, gainDB(high, 1), gainSlackDB, "norm %d, order %d", norm, order)
		}
	}
}

// A Bessel lowpass has no ripple: its gain falls monotonically from DC.
func TestBesselGainIsMonotonic(t *testing.T) {
	for _, norm := range besselNorms {
		for _, order := range []int{2, 5, 10} {
			f, err := Bessel(order, Lowpass(0.3), BesselNormalised(norm))
			require.NoError(t, err)
			previous := gainDB(f, 0)
			for i := 1; i < 500; i++ {
				gain := gainDB(f, float64(i)/500)
				assert.Less(t, gain, previous+gainSlackDB, "norm %d, order %d, w %g", norm, order, float64(i)/500)
				previous = gain
			}
		}
	}
}

// BesselNormDelay's prototype has a group delay of 1 s at DC. Moved to the
// prewarped edge W = 4·tan(π·wn/2) it's 1/W s, and the bilinear transform
// stretches frequency near DC by dΩ/dω = 2 at fs = 2, so the digital filter's
// delay at DC is exactly 2/W = 1/(2·tan(π·wn/2)) samples.
func TestBesselDelayAtDC(t *testing.T) {
	for _, wn := range []float64{0.01, 0.1, 0.3, 0.7} {
		for order := 1; order <= 10; order++ {
			f, err := Bessel(order, Lowpass(wn), BesselNormalised(BesselNormDelay))
			require.NoError(t, err)
			want := 1 / (2 * math.Tan(math.Pi*wn/2))
			assert.InEpsilon(t, want, groupDelayAtDC(f), 1e-12, "wn %g, order %d", wn, order)
		}
	}
}

// The prototype's group delay is maximally flat: it departs from its DC value
// as ω^(2N), so halving ω divides the departure by 2^(2N), within 5% from
// ω = 0.4. Butterworth's departs as ω², divided by about 4. Past order 5 the
// departure at ω = 0.2 is below rounding.
func TestBesselGroupDelayIsMaximallyFlat(t *testing.T) {
	// the analog group delay, Σ Re(1/(jω − p))
	delay := func(proto analogZPK, w float64) (d float64) {
		for _, p := range proto.Poles {
			d += real(1 / (complex(0, w) - p))
		}
		return
	}
	departure := func(proto analogZPK, w float64) float64 {
		return math.Abs(delay(proto, w)/delay(proto, 0) - 1)
	}
	const w = 0.4
	for order := 2; order <= 5; order++ {
		bessel := besselPrototype(order, BesselNormDelay)
		want := math.Pow(2, float64(2*order))
		assert.InEpsilon(t, want, departure(bessel, w)/departure(bessel, w/2), 0.05, "order %d", order)

		butter := butterworthPrototype(order)
		assert.Less(t, departure(butter, w)/departure(butter, w/2), 6.0, "order %d", order)
	}
}

// BesselNormPhase's prototype falls as 1/ω^N, as Butterworth's does: its
// poles' magnitudes multiply to 1, and its gain is 1.
func TestBesselPhasePrototypeMatchesButterworthAsymptote(t *testing.T) {
	for order := 1; order <= 12; order++ {
		proto := besselPrototype(order, BesselNormPhase)
		product := 1.0
		for _, p := range proto.Poles {
			product *= cmplx.Abs(p)
		}
		assert.InEpsilon(t, 1, product, 1e-14, "order %d", order)
		assert.InEpsilon(t, 1, proto.Gain, 1e-14, "order %d", order)
	}
}

// Order 1's one pole is at −1 with the -3 dB point at 1 rad/s, as
// Butterworth's, so BesselNormMagnitude gives exactly Butter's filter.
func TestBesselOrderOneIsButter(t *testing.T) {
	bessel, err := Bessel(1, Lowpass(0.3))
	require.NoError(t, err)
	butter, err := Butter(1, Lowpass(0.3))
	require.NoError(t, err)
	assert.Equal(t, butter, bessel)
}

func TestBesselPolesAreStable(t *testing.T) {
	bands := []Band{Lowpass(0.001), Lowpass(0.99), Highpass(0.0008), Bandpass(0.3, 0.31), Bandstop(0.01, 0.99)}
	for _, norm := range besselNorms {
		for _, band := range bands {
			for _, order := range []int{1, 2, 3, 8, 25} {
				f, err := Bessel(order, band, BesselNormalised(norm))
				require.NoError(t, err)
				for _, p := range f.Poles {
					assert.Less(t, cmplx.Abs(p), 1.0, "norm %d, band %v, order %d", norm, band, order)
				}
			}
		}
	}
}

func TestBesselRejectsInvalidInput(t *testing.T) {
	for _, order := range []int{0, -1} {
		_, err := Bessel(order, Lowpass(0.3))
		assert.ErrorIs(t, err, ErrInvalidOrder, "order %d", order)
	}
	_, err := Bessel(2, Band{})
	assert.ErrorIs(t, err, ErrInvalidBand)
	for _, norm := range []BesselNorm{-1, 3} {
		_, err := Bessel(2, Lowpass(0.3), BesselNormalised(norm))
		assert.ErrorIs(t, err, ErrInvalidMode, "BesselNorm %d", norm)
	}
}
