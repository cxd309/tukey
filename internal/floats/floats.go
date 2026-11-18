package floats

import (
	"math"
	"slices"
)

// Sum is a running sum with Neumaier's compensation:
// it keeps the low bits each addition round off, and adds them back in Total.
// A sum of many small values with a few large ones stays accurate,
// where a plain sum loses the small values' low bits.
// The zero Sum is empty.
type Sum struct {
	sum  float64
	lost float64 // the low bits rounded off so far
}

// Add adds v to the sum.
func (s *Sum) Add(v float64) {
	total := s.sum + v
	if math.Abs(s.sum) >= math.Abs(v) {
		s.lost += (s.sum - total) + v
	} else {
		s.lost += (v - total) + s.sum
	}
	s.sum = total
}

// Total is the sum of everything added, with the lost bits restored. An
// infinite or NaN sum is returned as it is, as its lost bits are then
// Inf - Inf, NaN, and mean nothing beside it.
func (s Sum) Total() (total float64) {
	if math.IsInf(s.sum, 0) || math.IsNaN(s.sum) {
		total = s.sum
		return
	}
	total = s.sum + s.lost
	return
}

// Scaled returns a new slice holding v multiplied by s
func Scaled(v []float64, s float64) (out []float64) {
	out = make([]float64, len(v))
	for i := range v {
		out[i] = v[i] * s
	}
	return
}

// Reversed returns a new slice holding v in reverse order
func Reversed(v []float64) (r []float64) {
	r = slices.Clone(v)
	slices.Reverse(r)
	return
}

// Dot is Σ w[k]·x[k], as if summed exactly and rounded once.
// Each product's rounding error is recovered by FMA, and every term is summed with compensation.
// The one exception is terms that cancel to far below their own size.
// An infinite product has no rounding error to recover, and is summed as it is.
func Dot(w, x []float64) (dot float64) {
	var sum Sum
	for k := range w {
		product := float64(w[k] * x[k])
		sum.Add(product)
		if !math.IsInf(product, 0) {
			sum.Add(math.FMA(w[k], x[k], -product))
		}
	}
	dot = sum.Total()
	return
}

// Midpoint is (a + b)/2, rounded once.
// a + b rounds once, and halving it is exact unless the half is subnormal.
// A subnormal half means a + b was on the subnormals' grid, so it was exact.
// Only an overflowing a + b needs a/2 + b/2, where each half is exact.
func Midpoint(a, b float64) (mid float64) {
	sum := a + b
	if math.IsInf(sum, 0) && !math.IsInf(a, 0) && !math.IsInf(b, 0) {
		mid = a/2 + b/2
		return
	}
	mid = sum / 2
	return
}

// InfiniteEnd is the value fraction t of the way along the line from a to b,
// when a or b is infinite, as its limit.
// t may be outside [0, 1], to extrapolate.
// With one end infinite, the line is that infinity on its side of the finite end,
// the finite end itself, and the opposite infinity beyond it.
// With both ends infinite it's a at t = 0, b at t = 1, and the shared infinity between them,
// and NaN anywhere else.
// Plain interpolation gives NaN instead, from Inf - Inf or 0·Inf.
// ok is false when neither end is infinite, for the caller to interpolate as usual.
func InfiniteEnd(a, b, t float64) (v float64, ok bool) {
	aInf, bInf := math.IsInf(a, 0), math.IsInf(b, 0)
	switch {
	case !aInf && !bInf:
		return 0, false
	case t == 0:
		return a, true
	case t == 1:
		return b, true
	case aInf && bInf:
		if a == b && t > 0 && t < 1 {
			return a, true
		}
		return math.NaN(), true
	case aInf && t < 1:
		return a, true
	case aInf:
		return -a, true
	case t > 0:
		return b, true
	}
	return -b, true
}
