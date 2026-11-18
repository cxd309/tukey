package special

import "math"

// SinCosPi is sin(π·x) and cos(π·x), within a few ULP for any x.
// x is reduced exactly to a quarter turn either side of a multiple of π/2,
// so an x far from 0 is as accurate as one near it.
func SinCosPi(x float64) (sin, cos float64) {
	n := math.Round(2 * x)
	y := x - n/2 // exact, and |y| ≤ 1/4
	s, c := math.Sincos(math.Pi * y)
	switch int64(n) & 3 {
	case 0:
		return s, c
	case 1:
		return c, -s
	case 2:
		return -s, -c
	default:
		return -c, s
	}
}

// SinCosPiProduct is sin(π·a·b) and cos(π·a·b).
// The rounding error of a·b is added back by its first-order term,
// so a product in the hundreds, such as a long FIR filter's outer taps, keeps its digits.
func SinCosPiProduct(a, b float64) (sin, cos float64) {
	t := float64(a * b)     // rounded and never fused, so e is its error
	e := math.FMA(a, b, -t) // a·b - t, exactly
	sin, cos = SinCosPi(t)
	return sin + math.Pi*e*cos, cos - math.Pi*e*sin
}
