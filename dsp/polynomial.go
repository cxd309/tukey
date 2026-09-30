package dsp

import (
	"cmp"
	"fmt"
	"math"
	"math/cmplx"
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

// conjugateTolerance is SciPy's tolerance for treating a root as real, or two roots as conjugates
// 100 machine episolons, relative to the root's magnitude
const conjugateTolerance = 100 * 0x1p-52

// collapseConjugates replaces each conjugate pair with on root
// the member with positive imaginary part, averaged with the other's conjugate to cancel rounding
// pairs come first, sorted by real part then imaginary part,
// followed by the real roots sorted by value
// returns ErrInvalidZPK if a complex root has no conjugate
// like SciPy's _cplxreal
func collapseConjugates(roots []complex128) (collapsed []complex128, err error) {
	reals, upper, lower := splitByImag(sortedByRealThenImag(roots))
	if len(upper) != len(lower) {
		return nil, fmt.Errorf("%w: %d complex roots above the real axis but %d below; they must come in conjugate pairs",
			ErrInvalidZPK, len(upper), len(lower))
	}
	alignConjugates(upper, lower)

	for i := range upper {
		if cmplx.Abs(upper[i]-cmplx.Conj(lower[i])) > conjugateTolerance*cmplx.Abs(lower[i]) {
			return nil, fmt.Errorf("%w: %v has no matching conjugate", ErrInvalidZPK, upper[i])
		}
		collapsed = append(collapsed, (upper[i]+cmplx.Conj(lower[i]))/2)
	}
	collapsed = append(collapsed, reals...)
	return
}

// sortedByRealThenImag returns a copy of roots sorted by real part then by the magnitude of the imaginary part
func sortedByRealThenImag(roots []complex128) (sorted []complex128) {
	sorted = slices.Clone(roots)
	slices.SortStableFunc(sorted, func(x, y complex128) int {
		return cmp.Or(cmp.Compare(real(x), real(y)), compareImagMagnitude(x, y))
	})
	return
}

// splitByImag separates roots into real ones (stored with exactly zero imaginary part),
// those above the real axis, and those below, keeping their order
func splitByImag(roots []complex128) (reals, upper, lower []complex128) {
	for _, r := range roots {
		switch {
		case math.Abs(imag(r)) <= conjugateTolerance*cmplx.Abs(r):
			reals = append(reals, complex(real(r), 0))
		case imag(r) > 0:
			upper = append(upper, r)
		default:
			lower = append(lower, r)
		}
	}
	return
}

// alignConjugates reorders both halves by imaginary part within each run of
// (nearly) equal real part, so upper[i] lines up with the conjugate of lower[i]
// runs are found in upper and applied to both, as SciPy does
func alignConjugates(upper, lower []complex128) {
	start := 0
	for i := 1; i <= len(upper); i++ {
		sameRun := i < len(upper) &&
			real(upper[i])-real(upper[i-1]) <= conjugateTolerance*cmplx.Abs(upper[i-1])
		if sameRun {
			continue
		}
		slices.SortStableFunc(upper[start:i], compareImagMagnitude)
		slices.SortStableFunc(lower[start:i], compareImagMagnitude)
		start = i
	}
}

// compareImagMagnitude orders roots by the magnitude of their imaginary part
func compareImagMagnitude(x, y complex128) (order int) {
	order = cmp.Compare(math.Abs(imag(x)), math.Abs(imag(y)))
	return
}
