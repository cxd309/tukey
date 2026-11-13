package dsp

import (
	"fmt"
	"math"
	"math/big"

	"github.com/cxd309/tukey/internal/floats"
)

// delayConditionLimit is the condition number of B(e^(-jω)),
// the sum of its terms' magnitudes over its own, above which it's evaluated in big.Float.
// Below it float64 is within 1e-14 of the exact delay on every reference vector.
// Above it, float64 loses up to 0.4% on butter(8, 0.01)'s (b, a).
const delayConditionLimit = 1000

// delayPrec is the working precision of the big.Float evaluation, in bits:
// enough for a condition number of 1e60
const delayPrec = 256

// bigPi is π to 365 bits, beyond delayPrec
var bigPi, _ = new(big.Float).SetPrec(delayPrec).SetString(
	"3.14159265358979323846264338327950288419716939937510582097494459230781640628620899862803482534211706798214808651")

// delayGrid is freqz's grid, ω = π·i/n, at which a polynomial's delay is found
type delayGrid struct {
	n int
}

// angle is m·ω/π for m = (2k − N)/2 at ω = π·i/n,
// as an exact fraction of a half turn, num/(2n), with num reduced to [0, 4n)
func (g delayGrid) angle(k, degree, i int) (num int) {
	fourN := 4 * g.n
	num = ((2*k - degree) * i) % fourN
	if num < 0 {
		num += fourN
	}
	return
}

// polynomialDelay is the group delay in samples of Σ c_k e^(−jkω),
// the numerator or denominator of a filter, at ω = π·i/n.
// Leading zeros delay it by a whole sample each, and trailing zeros don't change it.
// With m_k = k − N/2 for degree N,
//
//	B(e^(−jω)) = e^(−jNω/2)·F(ω),  F(ω) = Σ c_k e^(−j·m_k·ω) = R + jI
//
// so its delay is N/2 − Im(F′/F).
// F's real and imaginary parts pair each coefficient with its mirror, c_k ± c_(N−k):
// a symmetric or antisymmetric polynomial,
// such as a linear-phase FIR or a Butterworth numerator,
// has I or R exactly 0 and delay exactly N/2,
// even at a zero on the unit circle,
// where that's the limit from either side.
func (g delayGrid) polynomialDelay(c []float64, i int) (delay float64) {
	for len(c) > 0 && c[0] == 0 {
		c = c[1:]
		delay++
	}
	for len(c) > 0 && c[len(c)-1] == 0 {
		c = c[:len(c)-1]
	}
	degree := len(c) - 1
	half := float64(degree) / 2
	if degree <= 0 {
		return
	}

	var re, im, dRe, dIm floats.Sum
	var magnitude float64
	symmetric, antisymmetric := true, true
	for k := 0; 2*k <= degree; k++ {
		sum, diff := c[k]+c[degree-k], c[k]-c[degree-k]
		if 2*k == degree {
			sum, diff = c[k], 0
		}
		symmetric = symmetric && diff == 0
		antisymmetric = antisymmetric && sum == 0
		sin, cos := sinCosPi(float64(g.angle(k, degree, i)) / float64(2*g.n))
		m := float64(2*k-degree) / 2
		re.Add(float64(sum * cos))
		im.Add(float64(-diff * sin))
		dRe.Add(float64(-sum * m * sin))
		dIm.Add(float64(-diff * m * cos))
		magnitude += math.Abs(sum) + math.Abs(diff)
	}
	if symmetric || antisymmetric {
		return delay + half
	}
	r, x := re.Total(), im.Total()
	if !(magnitude < delayConditionLimit*math.Hypot(r, x)) {
		return delay + g.bigPolynomialDelay(c, i)
	}
	return delay + half - (r*dIm.Total()-x*dRe.Total())/(r*r+x*x)
}

// bigPolynomialDelay is polynomialDelay without its leading and trailing
// zeros, in big.Float, for a point where float64 has too few digits: near a
// zero of F, or at one, where the delay is its limit from either side.
// There, if F's first nonzero derivative is the m-th, F^(m),
//
//	F ≈ F^(m)·h^m/m! + F^(m+1)·h^(m+1)/(m+1)!, so Im(F′/F) → Im(F^(m+1)/((m+1)·F^(m)))
//
// with every F^(j) = Σ c_k (−j·m_k)^j e^(−j·m_k·ω).
func (g delayGrid) bigPolynomialDelay(c []float64, i int) (delay float64) {
	degree := len(c) - 1
	sins, coss := make([]*big.Float, len(c)), make([]*big.Float, len(c))
	for k := range c {
		sins[k], coss[k] = g.bigSinCosPi(g.angle(k, degree, i))
	}
	// derivative is F^(j)'s real and imaginary parts
	derivative := func(j int) (re, im *big.Float) {
		re, im = newBig(0), newBig(0)
		for k, coeff := range c {
			// c_k·m_k^j·(−j)^j·(cos − j·sin)
			term := newBig(coeff)
			m := newBig(float64(2*k-degree) / 2)
			for range j {
				term.Mul(term, m)
			}
			tc, ts := newBig(0).Mul(term, coss[k]), newBig(0).Mul(term, sins[k])
			switch j % 4 { // (−j)^j = 1, −j, −1, j
			case 0:
				re.Add(re, tc)
				im.Sub(im, ts)
			case 1:
				re.Sub(re, ts)
				im.Sub(im, tc)
			case 2:
				re.Sub(re, tc)
				im.Add(im, ts)
			default:
				re.Add(re, ts)
				im.Add(im, tc)
			}
		}
		return
	}
	// a coefficient's size, below which F^(j) is zero, rounding apart
	var largest float64
	for _, coeff := range c {
		largest = max(largest, math.Abs(coeff))
	}
	zero := newBig(largest)
	zero.SetMantExp(zero, -(delayPrec - 32))

	fRe, fIm := derivative(0)
	order := 0
	for isZero(fRe, fIm, zero) && order < degree {
		order++
		fRe, fIm = derivative(order)
	}
	nextRe, nextIm := derivative(order + 1)
	// Im(next/(f·(order+1))) = (f_re·next_im − f_im·next_re) / ((order+1)·|f|²)
	num := newBig(0).Mul(fRe, nextIm)
	num.Sub(num, newBig(0).Mul(fIm, nextRe))
	den := newBig(0).Mul(fRe, fRe)
	den.Add(den, newBig(0).Mul(fIm, fIm))
	den.Mul(den, newBig(float64(order+1)))
	ratio, _ := num.Quo(num, den).Float64()
	return float64(degree)/2 - ratio
}

// bigSinCosPi is sin and cos of π·num/(2n), num in [0, 4n), to delayPrec bits:
// reduced in integers to a quarter turn either side of a multiple of π/2, then by Taylor series
func (g delayGrid) bigSinCosPi(num int) (sin, cos *big.Float) {
	quadrant := (2*num + g.n) / (2 * g.n) // num/n, rounded
	y := num - quadrant*g.n               // |y| ≤ n/2: |θ| ≤ π/4
	theta := newBig(float64(y))
	theta.Mul(theta, bigPi)
	theta.Quo(theta, newBig(float64(2*g.n)))

	sin, cos = newBig(0), newBig(0)
	term := newBig(1) // θ^k/k!
	limit := newBig(1)
	limit.SetMantExp(limit, -delayPrec-8)
	for k := 0; ; k++ {
		switch k % 4 {
		case 0:
			cos.Add(cos, term)
		case 1:
			sin.Add(sin, term)
		case 2:
			cos.Sub(cos, term)
		default:
			sin.Sub(sin, term)
		}
		term.Mul(term, theta)
		term.Quo(term, newBig(float64(k+1)))
		if newBig(0).Abs(term).Cmp(limit) < 0 {
			break
		}
	}
	switch quadrant % 4 {
	case 1:
		sin, cos = cos, sin.Neg(sin)
	case 2:
		sin, cos = sin.Neg(sin), cos.Neg(cos)
	case 3:
		sin, cos = cos.Neg(cos), sin
	}
	return
}

// newBig is x at delayPrec
func newBig(x float64) *big.Float {
	return new(big.Float).SetPrec(delayPrec).SetFloat64(x)
}

// isZero reports whether re + j·im is below zero in both parts
func isZero(re, im, zero *big.Float) bool {
	abs := func(x *big.Float) *big.Float { return newBig(0).Abs(x) }
	return abs(re).Cmp(zero) < 0 && abs(im).Cmp(zero) < 0
}

// GroupDelay returns the group delay of the filter b/a in samples,
// at n frequencies w evenly spaced from 0 to just below π, as Freqz's:
// how many samples each frequency's envelope is delayed by, −d(phase)/dω.
//
//	w, delay, err := dsp.GroupDelay(b, a, 512)
//
// It's computed from the coefficients exactly, to rounding,
// not by differentiating the phase.
// A symmetric or antisymmetric b or a, such as a linear-phase FIR or a Butterworth numerator,
// contributes exactly half its order.
// Where b or a is ill-conditioned, its poles or zeros crowded against the unit circle,
// the point is evaluated in big.Float.
// At a zero on the unit circle, where the phase jumps by π, the delay is its limit from either side.
//
// A high-order filter's (b, a) is ill-conditioned,
// SOSGroupDelay on its sections is faster and gives the same answer.
//
// Returns:
//   - ErrInvalidCoefficients for coefficients Freqz rejects
//   - ErrInvalidLength if n < 0
//
// Equivalent to:
//   - scipy.signal.group_delay((b, a), n), which is 0 at a zero on the unit
//     circle and loses digits near one: 3% off for butter(8, 0.01)
//   - matlab [delay, w] = grpdelay(b, a, n), which also gives the limit
func GroupDelay(b, a []float64, n int) (w, delay []float64, err error) {
	if err = validateCoefficients(b, a); err != nil {
		return nil, nil, err
	}
	if n < 0 {
		return nil, nil, fmt.Errorf("%w: n must be >= 0, got %d", ErrInvalidLength, n)
	}
	w = freqzPoints(n)
	delay = make([]float64, n)
	g := delayGrid{n}
	for i := range delay {
		delay[i] = g.polynomialDelay(b, i) - g.polynomialDelay(a, i)
	}
	return
}

// SOSGroupDelay is GroupDelay for a filter given as second-order sections:
// the sum of each section's delay, which is better conditioned than the product's.
//
// Returns:
//   - ErrInvalidSOS for sections SOSFreqz rejects
//   - ErrInvalidLength if n < 0
//
// Equivalent to:
//   - scipy.signal.group_delay of the sections' product, which SciPy has no
//     sections form of
//   - matlab [delay, w] = grpdelay(sos, n)
func SOSGroupDelay(sos SOS, n int) (w, delay []float64, err error) {
	if err = sos.validate(); err != nil {
		return nil, nil, err
	}
	if n < 0 {
		return nil, nil, fmt.Errorf("%w: n must be >= 0, got %d", ErrInvalidLength, n)
	}
	w = freqzPoints(n)
	delay = make([]float64, n)
	g := delayGrid{n}
	for i := range delay {
		for _, section := range sos {
			delay[i] += g.polynomialDelay(section[:3], i) - g.polynomialDelay(section[3:], i)
		}
	}
	return
}
