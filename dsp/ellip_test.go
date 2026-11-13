package dsp

import (
	"math"
	"math/cmplx"
	"testing"

	"github.com/cxd309/tukey/internal/special"
	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

// ellipSpecs are (ripple, attenuation) pairs from gentle to sharp. 0.01 dB
// against 20 dB puts v0 past half its quarter period, through
// jacobiInQuarter's complementary branch, and 6 dB against 10 dB has the
// smallest 1 − m of the conformance cases.
var ellipSpecs = []struct{ ripple, attenuation float64 }{
	{1, 40},
	{0.5, 80},
	{0.01, 20},
	{6, 10},
}

// An Ellip band edge is where the gain first drops below −ripple, as for Cheby1.
func TestEllipGainAtBandEdges(t *testing.T) {
	for _, spec := range ellipSpecs {
		for _, c := range chebyBands {
			for order := 1; order <= 6; order++ {
				f, err := Ellip(order, spec.ripple, spec.attenuation, c.band)
				require.NoError(t, err)
				for _, edge := range c.edges {
					assert.InDelta(t, -spec.ripple, gainDB(f, edge), gainSlackDB, "%v, %s, order %d, edge %g", spec, c.name, order, edge)
				}
			}
		}
	}
}

// The passband's peak is 0 dB, so an even order's DC gain is −ripple, and an
// odd order's is 0 dB.
func TestEllipDCGain(t *testing.T) {
	for _, spec := range ellipSpecs {
		for order := 1; order <= 6; order++ {
			f, err := Ellip(order, spec.ripple, spec.attenuation, Lowpass(0.3))
			require.NoError(t, err)
			want := 0.0
			if order%2 == 0 {
				want = -spec.ripple
			}
			assert.InDelta(t, want, gainDB(f, 0), gainSlackDB, "%v, order %d", spec, order)
		}
	}
}

// The passband ripples between 0 dB and −ripple. Once the gain first reaches
// −attenuation, the stopband never rises above it, and its ripples' peaks
// touch it: equal ripple in both bands.
func TestEllipEqualRippleInBothBands(t *testing.T) {
	const cutoff, samples = 0.3, 20000
	for _, spec := range ellipSpecs {
		for _, order := range []int{2, 3, 4, 5} {
			f, err := Ellip(order, spec.ripple, spec.attenuation, Lowpass(cutoff))
			require.NoError(t, err)
			for i := range 500 {
				w := cutoff * float64(i) / 499
				gain := gainDB(f, w)
				assert.LessOrEqual(t, gain, gainSlackDB, "%v, order %d, w %g", spec, order, w)
				assert.GreaterOrEqual(t, gain, -spec.ripple-gainSlackDB, "%v, order %d, w %g", spec, order, w)
			}
			inStopband, peak := false, math.Inf(-1)
			for i := range samples {
				w := cutoff + (1-cutoff)*float64(i)/samples
				gain := gainDB(f, w)
				if !inStopband {
					inStopband = gain <= -spec.attenuation
					continue
				}
				assert.LessOrEqual(t, gain, -spec.attenuation+gainSlackDB, "%v, order %d, w %g", spec, order, w)
				peak = max(peak, gain)
			}
			require.True(t, inStopband, "%v, order %d never reaches the stopband", spec, order)
			// the sampled peaks fall just short of the true ones
			assert.InDelta(t, -spec.attenuation, peak, 1e-3, "%v, order %d", spec, order)
		}
	}
}

// The analog zeros are on the imaginary axis, so the digital ones are on the
// unit circle: true notches. An odd order's zero at infinity lands at Nyquist.
func TestEllipZerosAreOnUnitCircle(t *testing.T) {
	for _, order := range []int{2, 3, 6, 7} {
		f, err := Ellip(order, 1, 40, Lowpass(0.3))
		require.NoError(t, err)
		require.Len(t, f.Zeros, order)
		for _, z := range f.Zeros {
			assert.InDelta(t, 1, cmplx.Abs(z), 1e-14, "order %d, zero %v", order, z)
		}
		if order%2 == 1 {
			assert.Contains(t, f.Zeros, complex(-1, 0), "order %d", order)
		}
	}
}

// A first-order elliptic filter has no stopband ripple: it's Cheby1's.
func TestEllipOrderOneIsCheby1(t *testing.T) {
	for _, c := range chebyBands {
		ellip, err := Ellip(1, 1, 40, c.band)
		require.NoError(t, err)
		cheby, err := Cheby1(1, 1, c.band)
		require.NoError(t, err)
		assert.Equal(t, cheby, ellip, c.name)
	}
}

func TestEllipPolesAreStable(t *testing.T) {
	bands := []Band{Lowpass(0.001), Lowpass(0.99), Highpass(0.0008), Bandpass(0.3, 0.31), Bandstop(0.01, 0.99)}
	for _, band := range bands {
		for _, order := range []int{1, 2, 3, 8} {
			for _, spec := range ellipSpecs {
				f, err := Ellip(order, spec.ripple, spec.attenuation, band)
				require.NoError(t, err)
				for _, p := range f.Poles {
					assert.Less(t, cmplx.Abs(p), 1.0, "band %v, order %d, %v", band, order, spec)
				}
			}
		}
	}
}

// Either side of K/2, jacobiInQuarter agrees with Jacobi at a parameter
// where Jacobi is accurate across the whole quarter period.
func TestJacobiInQuarterMatchesJacobi(t *testing.T) {
	const mc = 0.3
	quarter := special.EllipKM1(mc)
	for _, fraction := range []float64{0.1, 0.49, 0.51, 0.9} {
		u := fraction * quarter
		sn, cn, dn := jacobiInQuarter(u, quarter-u, 1-mc, mc)
		wantSN, wantCN, wantDN := special.Jacobi(u, 1-mc, mc)
		assert.InEpsilon(t, wantSN, sn, 4e-15, "%gK", fraction)
		assert.InEpsilon(t, wantCN, cn, 4e-15, "%gK", fraction)
		assert.InEpsilon(t, wantDN, dn, 4e-15, "%gK", fraction)
	}
}

// Near K, at a parameter near 1, cn is small: Jacobi loses its relative
// digits there, but jacobiInQuarter keeps them. At u = K − δ, cn is
// √mc·sn(δ)/dn(δ) ≈ √mc·δ.
func TestJacobiInQuarterKeepsDigitsNearK(t *testing.T) {
	const mc, delta = 1e-6, 1e-3
	quarter := special.EllipKM1(mc)
	snDelta, _, dnDelta := special.Jacobi(delta, 1-mc, mc)
	want := math.Sqrt(mc) * snDelta / dnDelta
	_, cn, _ := jacobiInQuarter(quarter-delta, delta, 1-mc, mc)
	assert.InEpsilon(t, want, cn, 1e-15)
	_, direct, _ := special.Jacobi(quarter-delta, 1-mc, mc)
	assert.Greater(t, math.Abs(direct-want)/want, 1e-12, "Jacobi should drift here, or the case tests nothing")
}

func TestEllipRejectsInvalidInput(t *testing.T) {
	for _, dB := range []float64{0, -1, nan, math.Inf(1), math.Inf(-1)} {
		_, err := Ellip(2, dB, 40, Lowpass(0.3))
		assert.ErrorIs(t, err, ErrInvalidParameter, "rippleDB %g", dB)
		_, err = Ellip(2, 1, dB, Lowpass(0.3))
		assert.ErrorIs(t, err, ErrInvalidParameter, "attenuationDB %g", dB)
	}
	for _, attenuation := range []float64{1, 0.5} {
		_, err := Ellip(2, 1, attenuation, Lowpass(0.3))
		assert.ErrorIs(t, err, ErrInvalidParameter, "attenuationDB %g against 1 dB ripple", attenuation)
	}
	for _, order := range []int{0, -1} {
		_, err := Ellip(order, 1, 40, Lowpass(0.3))
		assert.ErrorIs(t, err, ErrInvalidOrder, "order %d", order)
	}
	_, err := Ellip(2, 1, 40, Band{})
	assert.ErrorIs(t, err, ErrInvalidBand)
}

// F(φ) + F(ψ) = K when tan φ · tan ψ = 1/√(1 − m), so the inverse sc and its
// complement, found separately, sum to the quarter period.
func TestIncompleteFComplementsSumToK(t *testing.T) {
	for _, spec := range ellipSpecs {
		epsilonSq := math.Expm1(spec.ripple * math.Ln10 / 10)
		stopSq := math.Expm1(spec.attenuation * math.Ln10 / 10)
		m1 := epsilonSq / stopSq
		sum := incompleteF(1/epsilonSq, m1) + incompleteF(stopSq, m1)
		assert.InEpsilon(t, special.EllipKM1(m1), sum, 1e-15, "%v", spec)
	}
}
