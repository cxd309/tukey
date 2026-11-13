package poly

import (
	"cmp"
	"fmt"
	"math"
	"math/big"
	"math/cmplx"
	"slices"
)

// ReverseBesselScale is how ReverseBesselRoots scales the roots
type ReverseBesselScale int

const (
	// UnitDelay leaves the roots as they are:
	// θ(0)/θ(s) has a group delay of 1 at DC
	UnitDelay ReverseBesselScale = iota
	// UnitProduct divides them by θ(0)^(1/N),
	// so their magnitudes multiply to 1
	UnitProduct
	// UnitHalfPower divides them by the frequency where |θ(0)/θ(jω)|² = 1/2,
	// moving it to 1
	UnitHalfPower
)

// besselStepLimit is how small Aberth's and Newton's last steps must be,
// relative to their root:
// far below float64's precision, so the result rounds same as the exact root
const besselStepLimit = 0x1p-80

// besselMaxIterations bounds the iterations.
// Past it the roots have failed to converge, which is a bug:
// with too little working precision for the order,
// rounding stops the steps shrinking below besselStepLimit,
// so it panics rather than return inaccurate roots.
const besselMaxIterations = 200

// ReverseBesselRoots returns the roots of the reverse Bessel polynomial θ of
// order N, scaled by scale, each rounded once from its exact value:
//
//	θ(s) = Σ (2N − k)! / (2^(N − k)·k!·(N − k)!) · s^k, for k = 0...N
//
// They're very sensitive to rounding:
// about 1.8 bits per order, so 1e13 at order 25, and any float64 evaluation of θ,
// by its coefficients or a recurrence, misses them by that factor times machine epsilon.
// So they're found in big.Float, at 2 bits per order beyond float64's 53 and a margin,
// by Aberth's method, with θ and θ′ from the recurrence
//
//	θ_n(s) = (2n − 1)·θ_(n−1)(s) + s²·θ_(n−2)(s),  θ_0 = 1, θ_1 = s + 1
//
// The roots come in conjugate pairs, plus one real root for an odd order,
// sorted by imaginary part. order < 1 panics.
//
// Equivalent to the poles of scipy.signal.besselap(N, norm) for norm
// "delay", "phase" and "mag", which evaluates θ with scipy.special.kve.
func ReverseBesselRoots(order int, scale ReverseBesselScale) (roots []complex128) {
	if order < 1 {
		panic(fmt.Sprintf("poly: reverse Bessel polynomial of order %d", order))
	}
	b := newBesselSolver(order)
	exact := b.roots()

	divisor := b.one()
	switch scale {
	case UnitProduct:
		divisor = b.nthRoot(b.theta0)
	case UnitHalfPower:
		divisor = b.halfPower()
	}

	roots = make([]complex128, order)
	for i, z := range exact {
		re, _ := new(big.Float).Quo(&z.re, divisor).Float64()
		im, _ := new(big.Float).Quo(&z.im, divisor).Float64()
		if math.Abs(im) < besselStepLimit*math.Abs(re) { // an odd order's real root
			im = 0
		}
		roots[i] = complex(re, im)
	}
	slices.SortFunc(roots, func(x, y complex128) int {
		return cmp.Or(cmp.Compare(imag(x), imag(y)), cmp.Compare(real(x), real(y)))
	})
	return
}

// besselSolver finds a reverse Bessel polynomial's roots at a working
// precision for its order
type besselSolver struct {
	order  int
	prec   uint
	theta0 *big.Float // θ(0) = (2N)! / (2^N·N!), exactly
}

func newBesselSolver(order int) (b besselSolver) {
	b.order = order
	b.prec = uint(53 + 64 + 2*order)
	theta0 := big.NewInt(1)
	for k := order + 1; k <= 2*order; k++ { // (2N)!/N!
		theta0.Mul(theta0, big.NewInt(int64(k)))
	}
	theta0.Rsh(theta0, uint(order))
	b.theta0 = new(big.Float).SetPrec(b.prec).SetInt(theta0)
	return
}

// one is 1 at the working precision
func (b besselSolver) one() (x *big.Float) {
	return new(big.Float).SetPrec(b.prec).SetInt64(1)
}

// theta is θ(s) and θ′(s), by the recurrence
func (b besselSolver) theta(s *bigComplex) (p, dp *bigComplex) {
	one := b.complex(1)
	prev, prevD := b.complex(1), b.complex(0)
	p, dp = b.complex(0).add(s, one), b.complex(1)
	sq := b.complex(0).mul(s, s)
	twoS := b.complex(0).add(s, s)
	t := b.complex(0)
	for n := 2; n <= b.order; n++ {
		m := new(big.Float).SetInt64(int64(2*n - 1))
		// θ′_n = m·θ′_(n−1) + 2s·θ_(n−2) + s²·θ′_(n−2)
		nextD := b.complex(0).mulReal(dp, m)
		nextD.add(nextD, t.mul(twoS, prev))
		nextD.add(nextD, t.mul(sq, prevD))
		// θ_n = m·θ_(n−1) + s²·θ_(n−2)
		next := b.complex(0).mulReal(p, m)
		next.add(next, t.mul(sq, prev))
		prev, prevD, p, dp = p, dp, next, nextD
	}
	return
}

// roots is θ's roots by Aberth's method, from a circle of radius θ(0)^(1/N),
// the roots' geometric mean, turned off the real axis
func (b besselSolver) roots() (z []*bigComplex) {
	n := float64(b.order)
	logTheta0 := lgamma(2*n+1) - n*math.Ln2 - lgamma(n+1)
	radius := math.Exp(logTheta0 / n)
	z = make([]*bigComplex, b.order)
	for k := range z {
		z[k] = b.complex(cmplx.Rect(radius, 2*math.Pi*float64(k)/n+0.4))
	}

	one := b.complex(1)
	sum, t, step := b.complex(0), b.complex(0), b.complex(0)
	for range besselMaxIterations {
		largest := 0.0
		for k := range z {
			p, dp := b.theta(z[k])
			ratio := b.complex(0).quo(p, dp)
			// step = ratio / (1 − ratio·Σ_(j≠k) 1/(z_k − z_j))
			sum.set(b.complex(0))
			for j := range z {
				if j != k {
					sum.add(sum, t.quo(one, t.sub(z[k], z[j])))
				}
			}
			step.quo(ratio, t.sub(one, t.mul(ratio, sum)))
			z[k].sub(z[k], step)
			largest = max(largest, cmplx.Abs(step.complex128())/cmplx.Abs(z[k].complex128()))
		}
		if largest < besselStepLimit {
			return
		}
	}
	panic(fmt.Sprintf("poly: reverse Bessel roots of order %d didn't converge", b.order))
}

// nthRoot is x^(1/N), by Newton's method from its float64 estimate
func (b besselSolver) nthRoot(x *big.Float) (r *big.Float) {
	estimate, _ := x.Float64()
	r = new(big.Float).SetPrec(b.prec).SetFloat64(math.Pow(estimate, 1/float64(b.order)))
	n := new(big.Float).SetPrec(b.prec).SetInt64(int64(b.order))
	for range besselMaxIterations {
		// r −= (r^N − x) / (N·r^(N−1))
		power := b.one()
		for range b.order - 1 {
			power.Mul(power, r)
		}
		residual := new(big.Float).SetPrec(b.prec).Mul(power, r)
		residual.Sub(residual, x)
		step := residual.Quo(residual, power.Mul(power, n))
		r.Sub(r, step)
		if ratio, _ := new(big.Float).Quo(step, r).Float64(); math.Abs(ratio) < besselStepLimit {
			return
		}
	}
	panic(fmt.Sprintf("poly: %d-th root didn't converge", b.order))
}

// halfPower is the ω where |θ(jω)|² = 2·θ(0)², so |θ(0)/θ(jω)|² = 1/2.
// |θ(jω)|² is a polynomial in ω² with positive coefficients,
// increasing and convex for ω > 0,
// so Newton's method from above it falls monotonically onto it.
func (b besselSolver) halfPower() (w *big.Float) {
	target := new(big.Float).SetPrec(b.prec).Mul(b.theta0, b.theta0)
	target.Add(target, target)
	w = new(big.Float).SetPrec(b.prec).SetFloat64(math.Sqrt(2*float64(b.order)) + 1)
	for range besselMaxIterations {
		s := b.complex(0)
		s.im.Set(w)
		p, dp := b.theta(s)
		// f = |θ|² − target, f′ = d/dω |θ(jω)|² = 2·Re(conj(θ)·j·θ′) = 2·(θi·θ′r − θr·θ′i)
		f := b.absSq(p)
		f.Sub(f, target)
		slope := new(big.Float).SetPrec(b.prec).Mul(&p.im, &dp.re)
		slope.Sub(slope, new(big.Float).SetPrec(b.prec).Mul(&p.re, &dp.im))
		slope.Add(slope, slope)
		step := f.Quo(f, slope)
		w.Sub(w, step)
		if ratio, _ := new(big.Float).Quo(step, w).Float64(); math.Abs(ratio) < besselStepLimit {
			return
		}
	}
	panic(fmt.Sprintf("poly: reverse Bessel half-power frequency of order %d didn't converge", b.order))
}

// absSq is |z|²
func (b besselSolver) absSq(z *bigComplex) (x *big.Float) {
	x = new(big.Float).SetPrec(b.prec).Mul(&z.re, &z.re)
	return x.Add(x, new(big.Float).SetPrec(b.prec).Mul(&z.im, &z.im))
}

// complex is c at the working precision
func (b besselSolver) complex(c complex128) (z *bigComplex) {
	z = &bigComplex{}
	z.re.SetPrec(b.prec).SetFloat64(real(c))
	z.im.SetPrec(b.prec).SetFloat64(imag(c))
	return
}

// lgamma is ln Γ(x) for x > 0
func lgamma(x float64) (y float64) {
	y, _ = math.Lgamma(x)
	return
}

// bigComplex is a complex number of big.Floats. Each method sets z to its
// result and returns it; z may be an argument, as for big.Float.
type bigComplex struct {
	re, im big.Float
}

func (z *bigComplex) set(x *bigComplex) *bigComplex {
	z.re.Set(&x.re)
	z.im.Set(&x.im)
	return z
}

func (z *bigComplex) add(x, y *bigComplex) *bigComplex {
	z.re.Add(&x.re, &y.re)
	z.im.Add(&x.im, &y.im)
	return z
}

func (z *bigComplex) sub(x, y *bigComplex) *bigComplex {
	z.re.Sub(&x.re, &y.re)
	z.im.Sub(&x.im, &y.im)
	return z
}

func (z *bigComplex) mul(x, y *bigComplex) *bigComplex {
	prec := max(z.re.Prec(), x.re.Prec())
	var rr, ii, ri, ir big.Float
	rr.SetPrec(prec).Mul(&x.re, &y.re)
	ii.SetPrec(prec).Mul(&x.im, &y.im)
	ri.SetPrec(prec).Mul(&x.re, &y.im)
	ir.SetPrec(prec).Mul(&x.im, &y.re)
	z.re.Sub(&rr, &ii)
	z.im.Add(&ri, &ir)
	return z
}

func (z *bigComplex) mulReal(x *bigComplex, y *big.Float) *bigComplex {
	z.re.Mul(&x.re, y)
	z.im.Mul(&x.im, y)
	return z
}

func (z *bigComplex) quo(x, y *bigComplex) *bigComplex {
	prec := max(z.re.Prec(), x.re.Prec())
	var denominator, t big.Float
	denominator.SetPrec(prec).Mul(&y.re, &y.re)
	denominator.Add(&denominator, t.SetPrec(prec).Mul(&y.im, &y.im))
	conj := bigComplex{}
	conj.re.SetPrec(prec).Set(&y.re)
	conj.im.SetPrec(prec).Neg(&y.im)
	z.mul(x, &conj)
	z.re.Quo(&z.re, &denominator)
	z.im.Quo(&z.im, &denominator)
	return z
}

func (z *bigComplex) complex128() complex128 {
	re, _ := z.re.Float64()
	im, _ := z.im.Float64()
	return complex(re, im)
}
