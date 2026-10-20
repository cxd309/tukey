// Package special is the special functions tukey needs
package special

import "math"

// Chebyshev coefficients of I0, highest order first and the constant term
// last, as chebyshev expects; from Cephes's i0, via SciPy's xsf
// (see THIRD_PARTY_NOTICES)
var (
	// exp(-x) I0(x) on [0, 8], in x/2 - 2; its limit as x -> 0 is 1
	i0Small = [...]float64{
		-4.41534164647933937950e-18, 3.33079451882223809783e-17, -2.43127984654795469359e-16,
		1.71539128555513303061e-15, -1.16853328779934516808e-14, 7.67618549860493561688e-14,
		-4.85644678311192946090e-13, 2.95505266312963983461e-12, -1.72682629144155570723e-11,
		9.67580903537323691224e-11, -5.18979560163526290666e-10, 2.65982372468238665035e-9,
		-1.30002500998624804212e-8, 6.04699502254191894932e-8, -2.67079385394061173391e-7,
		1.11738753912010371815e-6, -4.41673835845875056359e-6, 1.64484480707288970893e-5,
		-5.75419501008210370398e-5, 1.88502885095841655729e-4, -5.76375574538582365885e-4,
		1.63947561694133579842e-3, -4.32430999505057594430e-3, 1.05464603945949983183e-2,
		-2.37374148058994688156e-2, 4.93052842396707084878e-2, -9.49010970480476444210e-2,
		1.71620901522208775349e-1, -3.04682672343198398683e-1, 6.76795274409476084995e-1,
	}

	// exp(-x) sqrt(x) I0(x) on (8, inf), in 32/x - 2; its limit as x -> inf is 1/sqrt(2 pi)
	i0Large = [...]float64{
		-7.23318048787475395456e-18, -4.83050448594418207126e-18, 4.46562142029675999901e-17,
		3.46122286769746109310e-17, -2.82762398051658348494e-16, -3.42548561967721913462e-16,
		1.77256013305652638360e-15, 3.81168066935262242075e-15, -9.55484669882830764870e-15,
		-4.15056934728722208663e-14, 1.54008621752140982691e-14, 3.85277838274214270114e-13,
		7.18012445138366623367e-13, -1.79417853150680611778e-12, -1.32158118404477131188e-11,
		-3.14991652796324136454e-11, 1.18891471078464383424e-11, 4.94060238822496958910e-10,
		3.39623202570838634515e-9, 2.26666899049817806459e-8, 2.04891858946906374183e-7,
		2.89137052083475648297e-6, 6.88975834691682398426e-5, 3.36911647825569408990e-3,
		8.04490411014108831608e-1,
	}
)

// maxLog is the log of the largest float64: exp overflows above it
const maxLog = 7.097827128933839730962063185871e2

// I0 is the modified Bessel function of the first kind, order zero:
// the sum over k of ((x/2)^k / k!)^2
// it's even, I0(0) = 1, and it grows like exp(|x|) / sqrt(2 pi |x|),
// overflowing to +Inf beyond |x| ~ 713.99; I0(±Inf) is +Inf and I0(NaN) is NaN
//
// relative error within ~6e-16 (Cephes's measurement over [0, 30])
//
// the method and coefficients are Cephes's i0, via SciPy's xsf (see THIRD_PARTY_NOTICES):
// a Chebyshev expansion of exp(-x) I0(x) on [0, 8], and of exp(-x) sqrt(x) I0(x)
// in 32/x beyond; as current xsf, and unlike SciPy 1.18, it stays finite past
// exp's overflow at x ~ 709.78
func I0(x float64) (y float64) {
	x = math.Abs(x)
	switch {
	case x <= 8:
		y = math.Exp(x) * chebyshev(x/2-2, i0Small[:])
	case math.IsInf(x, 1):
		y = x // the formula below would give Inf/Inf
	case x > maxLog:
		// exp(x) overflows here though I0(x) is finite up to ~713.99,
		// so the exponential is applied in two halves
		e := math.Exp(x / 2)
		y = e * chebyshev(32/x-2, i0Large[:]) / math.Sqrt(x) * e
	default:
		y = math.Exp(x) * chebyshev(32/x-2, i0Large[:]) / math.Sqrt(x)
	}
	return
}

// chebyshev sums the Chebyshev series with coefficients c at x/2, by
// Clenshaw's recurrence; c is highest order first, constant term last,
// and the constant term is weighted by a half, as Cephes's chbevl
func chebyshev(x float64, c []float64) (y float64) {
	b0, b1, b2 := c[0], 0.0, 0.0
	for _, ck := range c[1:] {
		b2 = b1
		b1 = b0
		b0 = x*b1 - b2 + ck
	}
	y = 0.5 * (b0 - b2)
	return
}
