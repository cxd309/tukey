package special

import (
	"math"
	"testing"

	"github.com/stretchr/testify/assert"
)

// mpmath.ellipk at 60 digits
func TestEllipK(t *testing.T) {
	cases := []struct{ m, want float64 }{
		{0, 1.5707963267948966},
		{1e-10, 1.5707963268341665},
		{0.1, 1.6124413487202194},
		{0.5, 1.8540746773013719},
		{0.9, 2.5780921133481733},
		{0.999999, 8.294051463601063},
	}
	for _, c := range cases {
		assert.InEpsilon(t, c.want, EllipK(c.m), 4e-16, "K(%g)", c.m)
	}
}

// mpmath.ellipk(1 - p) at 60 digits (700 for 1e-300), where EllipK(1 - p)
// would round 1 - p: K(1 - 1e-20) would be K(1), +Inf
func TestEllipKM1(t *testing.T) {
	cases := []struct{ p, want float64 }{
		{1e-300, 346.77405831022674},
		{1e-20, 24.412145291060348},
		{1e-10, 12.899219826387599},
		{1e-3, 4.841132560550297},
	}
	for _, c := range cases {
		assert.InEpsilon(t, c.want, EllipKM1(c.p), 4e-16, "K(1 - %g)", c.p)
	}
	assert.Equal(t, math.Inf(1), EllipK(1))
	assert.Equal(t, math.Inf(1), EllipKM1(0))
	for _, m := range []float64{-0.1, 1.1, math.NaN()} {
		assert.True(t, math.IsNaN(EllipK(m)), "K(%g)", m)
	}
}

// mpmath.ellipfun at 60 digits, at parameter 1 - mc exactly
func TestJacobi(t *testing.T) {
	cases := []struct{ u, m, mc, sn, cn, dn float64 }{
		{0.3, 0.5, 0.5, 0.2934127331684554, 0.9559858618277871, 0.9782405041743613},
		{1.7, 0.9, 0.1, 0.9522305421080186, 0.30538008231819724, 0.4288721198784109},
		{2.5, 0.09999999999999998, 0.9, 0.6587664560557387, -0.7523474970887871, 0.9780606707344777},
		{10, 0.7, 0.3, 0.9780732803575124, 0.20826103392304346, 0.5747702678248819},
		{0.5, 9.999778782798785e-13, 0.999999999999, 0.47942553860418563, 0.8775825618903822, 0.9999999999998851},
		{0.4, 0.9999, 1e-4, 0.3799499045918845, 0.9250070648382195, 0.9250148680932744},
		// u = K/2, where cn = dn·sn = √(√mc / (1 + √mc)) and dn = mc^¼
		{2.995794670253498, 0.9999, 1e-4, 0.9950371902099892, 0.09950371902099893, 0.1},
	}
	for _, c := range cases {
		sn, cn, dn := Jacobi(c.u, c.m, c.mc)
		assert.InEpsilon(t, c.sn, sn, 2e-15, "sn(%g | %g)", c.u, c.m)
		assert.InEpsilon(t, c.cn, cn, 2e-15, "cn(%g | %g)", c.u, c.m)
		assert.InEpsilon(t, c.dn, dn, 2e-15, "dn(%g | %g)", c.u, c.m)
	}
}

// sn² + cn² = 1 and dn² + m·sn² = 1, needing no reference, up to 3 quarter periods
func TestJacobiIdentities(t *testing.T) {
	for _, mc := range []float64{0.9, 0.5, 0.1, 1e-3, 1e-6} {
		m := 1 - mc
		quarter := EllipKM1(mc)
		for _, fraction := range []float64{0.01, 0.1, 0.3, 0.5, 0.7, 0.9, 1, 1.5, 3} {
			u := fraction * quarter
			sn, cn, dn := Jacobi(u, m, mc)
			assert.InDelta(t, 1, sn*sn+cn*cn, 1e-15, "%gK, mc %g", fraction, mc)
			assert.InDelta(t, 1, dn*dn+m*sn*sn, 1e-15, "%gK, mc %g", fraction, mc)
		}
	}
}

// sn(K) = 1, cn(K) = 0 and dn(K) = √mc, at the quarter period
func TestJacobiAtQuarterPeriod(t *testing.T) {
	for _, mc := range []float64{0.9, 0.5, 0.1, 1e-3} {
		sn, cn, dn := Jacobi(EllipKM1(mc), 1-mc, mc)
		assert.InDelta(t, 1, sn, 1e-15, "mc %g", mc)
		assert.InDelta(t, 0, cn, 1e-15, "mc %g", mc)
		assert.InEpsilon(t, math.Sqrt(mc), dn, 1e-15, "mc %g", mc)
	}
}

// m = 0 and m = 1 are the circular and hyperbolic functions
func TestJacobiLimits(t *testing.T) {
	for _, u := range []float64{0.3, 1, 4} {
		sn, cn, dn := Jacobi(u, 0, 1)
		assert.Equal(t, [3]float64{math.Sin(u), math.Cos(u), 1}, [3]float64{sn, cn, dn})
		sn, cn, dn = Jacobi(u, 1, 0)
		assert.Equal(t, [3]float64{math.Tanh(u), 1 / math.Cosh(u), 1 / math.Cosh(u)}, [3]float64{sn, cn, dn})
	}
}

// mpmath.elliprf at 60 digits
func TestCarlsonRF(t *testing.T) {
	cases := []struct{ x, y, z, want float64 }{
		{1, 2, 0, 1.3110287771460598},
		{0.5, 1, 2, 0.9688576532724524},
		{1e-10, 1, 1, 1.5707863268734357},
		{2, 3, 4, 0.5840828416771517},
	}
	for _, c := range cases {
		assert.InEpsilon(t, c.want, CarlsonRF(c.x, c.y, c.z), 4e-16, "RF(%g, %g, %g)", c.x, c.y, c.z)
	}
}

// RF is symmetric, RF(x, x, x) = 1/√x, and K(m) = RF(0, 1 − m, 1): an
// independent check of EllipK by another method
func TestCarlsonRFIdentities(t *testing.T) {
	assert.InEpsilon(t, CarlsonRF(1, 2, 3), CarlsonRF(3, 1, 2), 4e-16)
	assert.InEpsilon(t, 1/math.Sqrt(7.0), CarlsonRF(7, 7, 7), 4e-16)
	for _, m := range []float64{0, 0.3, 0.9, 0.9999} {
		assert.InEpsilon(t, EllipK(m), CarlsonRF(0, 1-m, 1), 4e-16, "m %g", m)
	}
}

// mpmath.mfrom at 60 digits; e^−π is where m = ½
func TestParameterFromNome(t *testing.T) {
	cases := []struct{ q, want float64 }{
		{1e-6, 1.5999872000703995e-05},
		{0.01, 0.14787439154361495},
		{math.Exp(-math.Pi), 0.5000000000000001},
		{0.3, 0.9956043683608302},
	}
	for _, c := range cases {
		assert.InEpsilon(t, c.want, ParameterFromNome(c.q), 1e-15, "q %g", c.q)
	}
	assert.Equal(t, 0.0, ParameterFromNome(0))
	for _, q := range []float64{-0.1, 1, math.NaN()} {
		assert.True(t, math.IsNaN(ParameterFromNome(q)), "q %g", q)
	}
}

// above e^−π, m is found from the complementary nome; mpmath.mfrom at 60 digits
func TestParameterFromNomeNearOne(t *testing.T) {
	cases := []struct{ q, want float64 }{
		{0.5, 0.9999895221373104},
		{0.9, 1}, // 1 − m is 3e-40
		{0.999, 1},
	}
	for _, c := range cases {
		assert.InEpsilon(t, c.want, ParameterFromNome(c.q), 1e-15, "q %g", c.q)
	}
}

// either side of e^−π, the two branches meet at m = ½
func TestParameterFromNomeAtHalf(t *testing.T) {
	assert.InEpsilon(t, 0.5, ParameterFromNome(math.Nextafter(nomeAtHalf, 0)), 1e-15)
	assert.InEpsilon(t, 0.5000000000000001, ParameterFromNome(math.Nextafter(nomeAtHalf, 1)), 1e-15)
}
