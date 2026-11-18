package fft

import "math"

// twiddle returns e^(-2πi·k/n), the rotation every DFT is built from, accurately for any integer k.
//
// Computing it directly rounds the angle 2π·k/n before sin and cos see it,
// an absolute error that grows with the angle.
// Instead, k/n is reduced to an eighth of a turn by exact integer arithmetic,
// sin and cos are evaluated only on that small angle,
// and symmetry maps the results back.
// Quarter and half turns come out exact, as 1, -i, -1 and i.
func twiddle(k, n int) (w complex128) {
	k %= n
	if k < 0 {
		k += n
	}
	// the angle is (π/4)·(8k/n), so octant is which eighth of a turn it's in, and r how far through
	octant, r := (8*k)/n, (8*k)%n

	// the angle from the nearest multiple of π/4, below in an even octant and above in an odd one,
	// always in [0, π/4]
	frac := r
	if octant%2 == 1 {
		frac = n - r
	}
	s, c := math.Sincos(math.Pi / 4 * float64(frac) / float64(n))

	var cos, sin float64
	switch octant {
	case 0:
		cos, sin = c, s // φ
	case 1:
		cos, sin = s, c // π/2 - φ
	case 2:
		cos, sin = -s, c // π/2 + φ
	case 3:
		cos, sin = -c, s // π - φ
	case 4:
		cos, sin = -c, -s // π + φ
	case 5:
		cos, sin = -s, -c // 3π/2 - φ
	case 6:
		cos, sin = s, -c // 3π/2 + φ
	default:
		cos, sin = c, -s // 2π - φ
	}
	w = complex(cos, -sin) // e^(-iθ) = cos θ - i sin θ
	return
}

// twiddleTable returns e^(-2πi·j/n) for j = 0, 1, ..., n-1:
// every twiddle a transform of length n needs,
// or of any length dividing n, as w_m^j = w_n^(j·n/m).
func twiddleTable(n int) (table []complex128) {
	table = make([]complex128, n)
	for j := range table {
		table[j] = twiddle(j, n)
	}
	return
}
