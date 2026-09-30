package dsp

import (
	"fmt"
	"math"
	"math/cmplx"
	"slices"
)

// SOS is a filter as a cascade of second-order sections, applied in order
// each row is one section's coefficients b0, b1, b2, a0, a1, a2
// matched to SciPy's sos arrays
type SOS [][6]float64

// SOS converts the filter into second-order sections
// Each pole is paired with its nearest zero, starting from poles closes to unit circle
// these go in the last sections
// this keeps each section's peak gain small
// this makes SOS numerically robust at orders where BA's coefficients are not
//
// returns ErrInvalidZPK if complex zeros or poles aren't in conjugate pairs
//
// equivalent to scipy.signal.zpk2sos with its default pairing='nearest'
// the algorithm follows SciPy's implementation (see THIRD_PARTY_NOTICES)
func (f ZPK) SOS() (sos SOS, err error) {
	if len(f.Zeros) == 0 && len(f.Poles) == 0 {
		sos = SOS{{f.Gain, 0, 0, 1, 0, 0}}
		return
	}

	pr, sections, err := newPairing(f)
	if err != nil {
		return nil, err
	}

	sos = make(SOS, sections)
	// built last section first, so the poles closest to the unit circle end up last
	for si := sections - 1; si >= 0; si-- {
		zeros, poles := pr.next()
		if sos[si], err = section(zeros, poles); err != nil {
			return nil, err
		}
	}

	for i := range 3 {
		sos[0][i] *= f.Gain
	}
	return
}

// rootType selects which roots worstPole and nearestRoot consider
type rootType int

const (
	rootAny rootType = iota
	rootReal
	rootComplex
)

// matches reports whether r is of type t
func (t rootType) matches(r complex128) (ok bool) {
	switch t {
	case rootReal:
		ok = isReal(r)
	case rootComplex:
		ok = !isReal(r)
	default:
		ok = true
	}
	return
}

// isReal reports whether r has no imaginary part
// collapseConjugates stores real roots with an imaginary part of exactly zero
// matches numpy.isreal
func isReal(r complex128) (ok bool) {
	ok = imag(r) == 0
	return
}

// countReal is how many of roots are real
func countReal(roots []complex128) (n int) {
	for _, r := range roots {
		if isReal(r) {
			n++
		}
	}
	return
}

// worstPole is the index of the pole of type t closest to the unit circle
// the first on ties same as numpy.argmin
func worstPole(poles []complex128, t rootType) (index int) {
	index = -1
	best := math.Inf(1)
	for i, p := range poles {
		if d := math.Abs(1 - cmplx.Abs(p)); t.matches(p) && d < best {
			best, index = d, i
		}
	}
	if index < 0 {
		panic(fmt.Sprintf("dsp: no pole of type %d left to pair", t))
	}
	return
}

// nearestRoot is the index of the root of type t closest to target
// the first on ties, same as numpy's argsort
func nearestRoot(roots []complex128, target complex128, t rootType) (index int) {
	index = -1
	best := math.Inf(1)
	for i, r := range roots {
		if d := cmplx.Abs(r - target); t.matches(r) && d < best {
			best, index = d, i
		}
	}
	if index < 0 {
		panic(fmt.Sprintf("dsp: no zero of type %d left to pair", t))
	}
	return
}

// takeAt removes and returns roots[i], keeping the rest in order, same as numpy.delete
func takeAt(roots []complex128, i int) (r complex128, rest []complex128) {
	r = roots[i]
	rest = slices.Delete(roots, i, i+1)
	return
}

// padRoots returns a copy of roots extended to length n with roots at the origin
func padRoots(roots []complex128, n int) (padded []complex128) {
	padded = make([]complex128, n)
	copy(padded, roots)
	return
}

// section builds one second-order section from up to two zeros and two poles
// right-aligning each polynomial so a missing root leave a leading zero
// same as SciPy _single_zpksos
func section(zeros, poles []complex128) (s [6]float64, err error) {
	b, err := realCoeffs(polyFromRoots(zeros))
	if err != nil {
		return s, fmt.Errorf("%w: %v", ErrInvalidZPK, err)
	}
	a, err := realCoeffs(polyFromRoots(poles))
	if err != nil {
		return s, fmt.Errorf("%w: %v", ErrInvalidZPK, err)
	}
	copy(s[3-len(b):3], b)
	copy(s[6-len(a):6], a)
	return
}
