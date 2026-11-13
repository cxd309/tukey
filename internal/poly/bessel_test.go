package poly

import (
	"cmp"
	"math"
	"math/cmplx"
	"slices"
	"testing"

	"github.com/stretchr/testify/assert"
)

var besselScales = []ReverseBesselScale{UnitDelay, UnitProduct, UnitHalfPower}

// besselMaxTestOrder is the highest order the property tests run:
// higher orders are slow, and order 25 is checked against its exact roots
const besselMaxTestOrder = 20

// Orders 1 and 2 in closed form. θ_1 = s + 1 and θ_2 = s² + 3s + 3,
// whose roots are -3/2 ± j√3/2, of magnitude √3.
// Order 2's |θ(jω)|² = ω⁴ + 3ω² + 9 is 2·θ(0)² = 18 at ω² = 3(√5 − 1)/2.
func TestReverseBesselRootsLowOrders(t *testing.T) {
	assert.Equal(t, []complex128{-1}, ReverseBesselRoots(1, UnitDelay))
	assert.Equal(t, []complex128{-1}, ReverseBesselRoots(1, UnitProduct))
	assert.Equal(t, []complex128{-1}, ReverseBesselRoots(1, UnitHalfPower))

	half := math.Sqrt(3) / 2
	assert.Equal(t, []complex128{complex(-1.5, -half), complex(-1.5, half)}, ReverseBesselRoots(2, UnitDelay))
	assert.Equal(t, []complex128{complex(-half, -0.5), complex(-half, 0.5)}, ReverseBesselRoots(2, UnitProduct))

	cutoff := math.Sqrt(3 * (math.Sqrt(5) - 1) / 2)
	got := ReverseBesselRoots(2, UnitHalfPower)
	for i, want := range []complex128{complex(-1.5, -half), complex(-1.5, half)} {
		assert.InEpsilon(t, real(want)/cutoff, real(got[i]), 1e-15)
		assert.InEpsilon(t, imag(want)/cutoff, imag(got[i]), 1e-15)
	}
}

// Order 25, MATLAB's highest,
// where θ's roots are most sensitive to rounding of any order MATLAB tabulates:
// 8e12 times, so float64 Aberth on the coefficients misses them by 1e-4.
// These are mpmath's roots at 50 digits, scaled to UnitProduct, rounded once.
// They're matched exactly: big.Float arithmetic is the same on every platform.
func TestReverseBesselRootsAreExact(t *testing.T) {
	upper := []complex128{
		complex(-0.9062073871811709, 0.0),
		complex(-0.902883339022802, 0.09307713118510297),
		complex(-0.8928551459883549, 0.18630689698043007),
		complex(-0.8759497989677858, 0.2798521321771409),
		complex(-0.851861688655402, 0.3738977875907595),
		complex(-0.820122604393688, 0.46866685746569664),
		complex(-0.7800496278186497, 0.564444121034971),
		complex(-0.7306549271849968, 0.6616149647357749),
		complex(-0.6704827128029559, 0.760734885816784),
		complex(-0.5972898661335557, 0.8626676330388029),
		complex(-0.5073362861078469, 0.9689006305344868),
		complex(-0.393452987819108, 1.0824339271738317),
		complex(-0.2373280669322029, 1.2114766583825654),
	}
	var want []complex128
	for _, r := range slices.Backward(upper[1:]) {
		want = append(want, cmplx.Conj(r))
	}
	want = append(want, upper...)
	assert.Equal(t, want, ReverseBesselRoots(25, UnitProduct))
}

// Every order's roots are in the left half-plane, in exact conjugate pairs
// sorted by imaginary part, with one real root, exactly real, for an odd order.
func TestReverseBesselRootsArePairedAndStable(t *testing.T) {
	for _, scale := range besselScales {
		for order := 1; order <= besselMaxTestOrder; order++ {
			roots := ReverseBesselRoots(order, scale)
			assert.Len(t, roots, order)
			assert.True(t, slices.IsSortedFunc(roots, func(x, y complex128) int {
				return cmp.Compare(imag(x), imag(y))
			}), "scale %d, order %d", scale, order)
			for i, r := range roots {
				assert.Less(t, real(r), 0.0, "scale %d, order %d", scale, order)
				assert.Equal(t, cmplx.Conj(r), roots[order-1-i], "scale %d, order %d", scale, order)
			}
			if order%2 == 1 {
				assert.Zero(t, imag(roots[order/2]), "scale %d, order %d", scale, order)
			}
		}
	}
}

// The roots' sums by Vieta: θ's s^(N−1) coefficient is N(N + 1)/2 and its
// s¹ and s⁰ coefficients are equal, so Σ r = −N(N + 1)/2 and Σ 1/r = −1:
// the group delay at DC is 1.
func TestReverseBesselRootsUnitDelay(t *testing.T) {
	for order := 1; order <= besselMaxTestOrder; order++ {
		roots := ReverseBesselRoots(order, UnitDelay)
		var sum, reciprocals complex128
		for _, r := range roots {
			sum += r
			reciprocals += 1 / r
		}
		assert.InEpsilon(t, -float64(order*(order+1))/2, real(sum), 1e-14, "order %d", order)
		assert.InEpsilon(t, -1, real(reciprocals), 1e-14, "order %d", order)
	}
}

// UnitProduct's magnitudes multiply to 1.
func TestReverseBesselRootsUnitProduct(t *testing.T) {
	for order := 1; order <= besselMaxTestOrder; order++ {
		product := 1.0
		for _, r := range ReverseBesselRoots(order, UnitProduct) {
			product *= cmplx.Abs(r)
		}
		assert.InEpsilon(t, 1, product, 1e-14, "order %d", order)
	}
}

// UnitHalfPower's filter, Π(−r)/Π(s − r), has half its DC power at 1 rad/s.
func TestReverseBesselRootsUnitHalfPower(t *testing.T) {
	for order := 1; order <= besselMaxTestOrder; order++ {
		h := complex(1, 0)
		for _, r := range ReverseBesselRoots(order, UnitHalfPower) {
			h *= -r / (1i - r)
		}
		assert.InEpsilon(t, 0.5, real(h*cmplx.Conj(h)), 1e-14, "order %d", order)
	}
}

// Each scale divides UnitDelay's roots by one positive factor.
func TestReverseBesselScalesAreProportional(t *testing.T) {
	for order := 1; order <= besselMaxTestOrder; order++ {
		delay := ReverseBesselRoots(order, UnitDelay)
		for _, scale := range besselScales[1:] {
			scaled := ReverseBesselRoots(order, scale)
			factor := real(delay[0]) / real(scaled[0])
			for i := range delay {
				assert.InEpsilon(t, factor, real(delay[i])/real(scaled[i]), 1e-15, "scale %d, order %d", scale, order)
			}
		}
	}
}

func TestReverseBesselRootsPanicsBelowOrderOne(t *testing.T) {
	for _, order := range []int{0, -1} {
		assert.Panics(t, func() { ReverseBesselRoots(order, UnitDelay) }, "order %d", order)
	}
}
