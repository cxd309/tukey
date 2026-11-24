package special

import "math"

// DoubleDouble is the unevaluated sum Hi + Lo, with |Lo| at most half an ULP of Hi,
// which carries about 32 significant digits, twice a float64's.
// Its arithmetic is within a few parts in 2^104 of exact,
// for sums such as a chirp's phase, where float64 would leave too few digits after removing whole cycles.
//
// Algorithms from Hida, Li and Bailey, "Library for double-double and quad-double arithmetic" (2007),
// with products split by math.FMA.
type DoubleDouble struct {
	Hi, Lo float64
}

// Exact is x as a DoubleDouble.
func Exact(x float64) (d DoubleDouble) {
	d = DoubleDouble{Hi: x}
	return
}

// SumOf is a + b, exactly.
func SumOf(a, b float64) (d DoubleDouble) {
	hi := a + b
	v := hi - a
	d = DoubleDouble{Hi: hi, Lo: (a - (hi - v)) + (b - v)}
	return
}

// ProductOf is a·b, exactly unless it underflows.
func ProductOf(a, b float64) (d DoubleDouble) {
	hi := float64(a * b) // rounded and never fused, so the FMA below is its error
	d = DoubleDouble{Hi: hi, Lo: math.FMA(a, b, -hi)}
	return
}

// quickSum is a + b, exactly, for |a| ≥ |b| or a = 0.
func quickSum(a, b float64) (d DoubleDouble) {
	hi := a + b
	d = DoubleDouble{Hi: hi, Lo: b - (hi - a)}
	return
}

// Float64 is d rounded to a float64.
func (d DoubleDouble) Float64() (x float64) {
	x = d.Hi + d.Lo
	return
}

// Neg is -d.
func (d DoubleDouble) Neg() (n DoubleDouble) {
	n = DoubleDouble{Hi: -d.Hi, Lo: -d.Lo}
	return
}

// Add is d + e.
func (d DoubleDouble) Add(e DoubleDouble) (sum DoubleDouble) {
	s := SumOf(d.Hi, e.Hi)
	t := SumOf(d.Lo, e.Lo)
	sum = quickSum(s.Hi, s.Lo+t.Hi)
	sum = quickSum(sum.Hi, sum.Lo+t.Lo)
	return
}

// Sub is d - e.
func (d DoubleDouble) Sub(e DoubleDouble) (diff DoubleDouble) {
	diff = d.Add(e.Neg())
	return
}

// Mul is d·e.
func (d DoubleDouble) Mul(e DoubleDouble) (product DoubleDouble) {
	p := ProductOf(d.Hi, e.Hi)
	cross := float64(d.Hi*e.Lo) + float64(d.Lo*e.Hi)
	product = quickSum(p.Hi, p.Lo+cross)
	return
}

// Div is d/e, three float64 quotients each taken from the remainder of the last.
func (d DoubleDouble) Div(e DoubleDouble) (quotient DoubleDouble) {
	q1 := d.Hi / e.Hi
	r := d.Sub(e.Mul(Exact(q1)))
	q2 := r.Hi / e.Hi
	r = r.Sub(e.Mul(Exact(q2)))
	q3 := r.Hi / e.Hi
	quotient = quickSum(q1, q2).Add(Exact(q3))
	return
}

// ldexp is d·2^k, exactly unless it overflows or underflows.
func (d DoubleDouble) ldexp(k int) (scaled DoubleDouble) {
	scaled = DoubleDouble{Hi: math.Ldexp(d.Hi, k), Lo: math.Ldexp(d.Lo, k)}
	return
}

// ln2 is ln(2) as a DoubleDouble, from mpmath.
var ln2 = DoubleDouble{Hi: 0.6931471805599453, Lo: 2.3190468138462996e-17}

// expSquarings is how many times ExpM1 halves its argument before its series, and squares after.
// Each halving shortens the series, and each squaring doubles its relative error,
// so 10 keeps the series to 9 terms and the error within 2^-94.
const expSquarings = 10

// Exp is e^d, within about 2^-94 relative of exact.
// It's +Inf above ln(MaxFloat64), and 0 far below 0.
func Exp(d DoubleDouble) (e DoubleDouble) {
	k, s, ok := expParts(d)
	if !ok {
		return Exact(math.Exp(d.Hi))
	}
	e = s.Add(Exact(1)).ldexp(k)
	return
}

// ExpM1 is e^d - 1, within about 2^-94 relative of exact,
// and as accurate for a d near 0, where e^d - 1 cancels.
// It's +Inf above ln(MaxFloat64), and -1 far below 0.
func ExpM1(d DoubleDouble) (e DoubleDouble) {
	k, s, ok := expParts(d)
	switch {
	case !ok:
		return Exact(math.Expm1(d.Hi))
	case k == 0:
		return s
	}
	e = s.Add(Exact(1)).ldexp(k).Sub(Exact(1))
	return
}

// expParts splits e^d into 2^k·(1 + s), with |s| < ½,
// or isn't ok for a NaN d, or one too large or too small for a finite, normal result.
func expParts(d DoubleDouble) (k int, s DoubleDouble, ok bool) {
	if math.IsNaN(d.Hi) || d.Hi > 709.78 || d.Hi < -708.39 {
		return
	}
	// d = k·ln(2) + r, with |r| ≤ ln(2)/2, then r is halved expSquarings times
	whole := math.Round(d.Hi / ln2.Hi)
	r := d.Sub(ln2.Mul(Exact(whole))).ldexp(-expSquarings)

	// e^r - 1 by its series, to below 2^-110 of the sum
	term := r
	for n := 2.0; math.Abs(term.Hi) > 0x1p-110*math.Abs(s.Hi); n++ {
		s = s.Add(term)
		term = term.Mul(r).Div(Exact(n))
	}
	// e^(2r) - 1 = 2(e^r - 1) + (e^r - 1)², which keeps a small result's digits
	for range expSquarings {
		s = s.ldexp(1).Add(s.Mul(s))
	}
	k, ok = int(whole), true
	return
}

// Log is ln(d), within about 2^-100 relative of exact.
// d = m·2^k with m between √½ and √2, and ln(m) is math.Log's estimate y,
// corrected by one Newton step, y + m·e^-y - 1.
// It's -Inf at 0, and NaN below 0.
func Log(d DoubleDouble) (y DoubleDouble) {
	switch {
	case d.Hi == 0:
		return Exact(math.Inf(-1))
	case d.Hi < 0 || math.IsNaN(d.Hi):
		return Exact(math.NaN())
	case math.IsInf(d.Hi, 1):
		return d
	}
	_, k := math.Frexp(d.Hi)
	if math.Abs(math.Ldexp(d.Hi, -k)) < math.Sqrt2/2 {
		k--
	}
	m := d.ldexp(-k)
	y = Exact(math.Log(m.Hi))
	// m·e^-y - 1 = (m - 1) + m·(e^-y - 1), both small when m is near 1
	step := m.Sub(Exact(1)).Add(m.Mul(ExpM1(y.Neg())))
	y = y.Add(step).Add(ln2.Mul(Exact(float64(k))))
	return
}
