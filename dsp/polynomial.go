package dsp

import (
	"fmt"
	"math"
	"slices"
)

// polyFromRoots returns the coefficients [x^n, x^(n-1), ..., x^0]
// of prod(x - r) for the given roots
// same convention as numpy.poly
func polyFromRoots(roots []complex128) (coeffs []complex128) {
	coeffs = []complex128{1}
	for _, r := range roots {
		next := make([]complex128, len(coeffs)+1)
		next[0] = coeffs[0]
		for i := 1; i < len(coeffs); i++ {
			next[i] = coeffs[i] - r*coeffs[i-1]
		}
		next[len(coeffs)] = -r * coeffs[len(coeffs)-1]
		coeffs = next
	}
	return
}

// realCoeffs drops the imaginary part of polynomial coefficients
// that should be real because their roots came in conjugate pairs
// it returns an error if the imaginary part is > floating-point rounding noise
func realCoeffs(complexCoeffs []complex128) (coeffs []float64, err error) {
	coeffs = make([]float64, len(complexCoeffs))
	for i, v := range complexCoeffs {
		if math.Abs(imag(v)) > 1e-9*math.Max(math.Abs(real(v)), 1) {
			return nil, fmt.Errorf("coefficient %d has imaginary part %v, complex roots must come in conjugate pairs", i, imag(v))
		}
		coeffs[i] = real(v)
	}
	return
}

// mapRoots applies fn to every root, returning a new slice
func mapRoots(roots []complex128, fn func(complex128) complex128) (mapped []complex128) {
	mapped = make([]complex128, len(roots))
	for i, r := range roots {
		mapped[i] = fn(r)
	}
	return
}

// prod returns the product of all values, or 1 for an empty slice
func prod(values []complex128) (p complex128) {
	p = 1
	for _, v := range values {
		p *= v
	}
	return
}

// appendRepeated returns roots followed by n copies of root, leaving roots unmodified
// n < 0 panics, which flags a filter with more zeros than poles
func appendRepeated(roots []complex128, root complex128, n int) (padded []complex128) {
	padded = slices.Concat(roots, slices.Repeat([]complex128{root}, n))
	return
}
