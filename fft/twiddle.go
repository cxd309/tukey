package fft

import "math"

// twiddle returns e^(-2 pi i k/n), the rotation every DFT is build from,
// accurately for any integer k
//
// Computing it directly as e^(-2 pi i k/n) rounds the angle 2 pi kn before sin and cos see it
// this is an absolute error that grows with the angle
//
// Instead, k/n is reduced to an eigth of a turn using exact integer arithmetic
// sin and cos are evaluated only on that small and symmetry maps the results back

// Quarter and half turns come out exact: 1, -i, -1, i
func twiddle(k, n int) (w complex128) {
	k %= n
	if k < 0 {
		k += n
	}
	// the angle is (pi / 4)(8 k/n) which eighth of a turn it's in and how far through
	octant, r := (8*k)/n, (8*k)%n

	// the angle measured from the nearest multiple of pi/4 below (even octant) or above (odd octant)
	// always in [0, pi/4]
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

// twiddleTable returns e^(-2 pi i j/n) for j=0...n-1
// every twiddle a transform of length n,
// or of any length dividing n, needs, since w_m^j = w_n^(j n/m)
func twiddleTable(n int) (table []complex128) {
	table = make([]complex128, n)
	for j := range table {
		table[j] = twiddle(j, n)
	}
	return
}
