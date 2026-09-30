package poly

import (
	"cmp"
	"fmt"
	"math"
	"math/cmplx"
	"slices"
)

// conjugateTolerance is SciPy's tolerance for treating a root as real, or two roots as conjugates
// 100 machine epsilons, relative to the root's magnitude
const conjugateTolerance = 100 * 0x1p-52

// CollapseConjugates replaces each conjugate pair with one root:
// the member with positive imaginary part, averaged with the other's conjugate to cancel rounding
// pairs come first, sorted by real part then imaginary part,
// followed by the real roots sorted by value
// returns ErrUnpairedConjugate if a complex root has no conjugate
// like SciPy's _cplxreal
func CollapseConjugates(roots []complex128) (collapsed []complex128, err error) {
	reals, upper, lower := splitByImag(sortedByRealThenImag(roots))
	if len(upper) != len(lower) {
		return nil, fmt.Errorf("%w: %d complex roots above the real axis but %d below",
			ErrUnpairedConjugate, len(upper), len(lower))
	}
	alignConjugates(upper, lower)

	for i := range upper {
		if cmplx.Abs(upper[i]-cmplx.Conj(lower[i])) > conjugateTolerance*cmplx.Abs(lower[i]) {
			return nil, fmt.Errorf("%w: %v has no matching conjugate", ErrUnpairedConjugate, upper[i])
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
