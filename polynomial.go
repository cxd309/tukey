package dsp

import (
	"fmt"
	"math"
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
// it panics if the imaginary part isn't just floating-point rounding noise
// since that would mean a bug earlier in the pipeline not a caller error
func realCoeffs(complexCoeffs []complex128) (coeffs []float64) {
	coeffs = make([]float64, len(complexCoeffs))
	for i, v := range complexCoeffs {
		if math.Abs(imag(v)) > 1e-9*math.Max(math.Abs(real(v)), 1) {
			panic(fmt.Sprintf("dsp: coefficient %d has non-negligible imaginary part: %v", i, v))
		}
		coeffs[i] = real(v)
	}
	return
}
