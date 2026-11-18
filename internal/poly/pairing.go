package poly

import (
	"fmt"
	"math"
	"math/cmplx"
	"slices"
)

// Pairing holds the zeros and poles not yet placed in a second-order section.
// Each conjugate pair is stored once, see CollapseConjugates.
// It follows the loop in SciPy's zpk2sos, with its default pairing='nearest'.
type Pairing struct {
	zeros, poles []complex128
}

// NewPairing balances zeros and poles, collapses their conjugate pairs,
// and returns the pool and how many sections it makes.
//
// Returns: ErrUnpairedConjugate, wrapped, if complex zeros or poles aren't in conjugate pairs.
func NewPairing(zeros, poles []complex128) (pr Pairing, sections int, err error) {
	// The counts are balanced with roots at the origin.
	// An odd count gets one more of each, so every section has two zeros and two poles.
	n := max(len(zeros), len(poles))
	zeros, poles = padRoots(zeros, n), padRoots(poles, n)
	sections = (n + 1) / 2
	if n%2 == 1 {
		zeros, poles = append(zeros, 0), append(poles, 0)
	}
	if pr.zeros, err = CollapseConjugates(zeros); err != nil {
		return pr, 0, fmt.Errorf("zeros: %w", err)
	}
	if pr.poles, err = CollapseConjugates(poles); err != nil {
		return pr, 0, fmt.Errorf("poles: %w", err)
	}
	return
}

// Next removes and returns the zeros and poles for the next section,
// starting from the pole closest to the unit circle.
func (pr *Pairing) Next() (zeros, poles []complex128) {
	p1 := pr.takePole(rootAny)
	switch {
	case isReal(p1) && countReal(pr.poles) == 0:
		zeros, poles = pr.lastRealPole(p1)
	case pr.mustKeepRealPairTogether(p1):
		zeros, poles = pr.complexPair(p1)
	default:
		zeros, poles = pr.secondOrder(p1)
	}
	return
}

// rootType selects which roots worstPole and nearestRoot consider.
type rootType int

const (
	rootAny rootType = iota
	rootReal
	rootComplex
)

// matches reports whether r is of type t.
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

// lastRealPole pairs the only remaining real pole with the nearest real zero,
// plus a pole and zero at the origin.
func (pr *Pairing) lastRealPole(p1 complex128) (zeros, poles []complex128) {
	zeros = []complex128{pr.takeZero(p1, rootReal), 0}
	poles = []complex128{p1, 0}
	return
}

// mustKeepRealPairTogether reports whether exactly one real pole and one real zero remain for a later section,
// in which case the complex pole p1 must take a complex zero.
func (pr *Pairing) mustKeepRealPairTogether(p1 complex128) (must bool) {
	must = !isReal(p1) && len(pr.poles)+1 == len(pr.zeros) &&
		countReal(pr.poles) == 1 && countReal(pr.zeros) == 1
	return
}

// complexPair pairs the complex pole p1 with the nearest complex zero, and both with their conjugates.
func (pr *Pairing) complexPair(p1 complex128) (zeros, poles []complex128) {
	z1 := pr.takeZero(p1, rootComplex)
	zeros = []complex128{z1, cmplx.Conj(z1)}
	poles = []complex128{p1, cmplx.Conj(p1)}
	return
}

// secondOrder completes a section around p1:
// its conjugate, or if it's real, the next real pole closest to the unit circle,
// then the nearest zero and whatever completes it,
// which is its conjugate if it's complex, or the nearest other real zero.
func (pr *Pairing) secondOrder(p1 complex128) (zeros, poles []complex128) {
	poles = []complex128{p1, cmplx.Conj(p1)}
	if isReal(p1) {
		poles[1] = pr.takePole(rootReal)
	}
	if len(pr.zeros) == 0 {
		return
	}

	z1 := pr.takeZero(p1, rootAny)
	switch {
	case !isReal(z1):
		zeros = []complex128{z1, cmplx.Conj(z1)}
	case len(pr.zeros) > 0:
		zeros = []complex128{z1, pr.takeZero(p1, rootReal)}
	default:
		zeros = []complex128{z1}
	}
	return
}

// takePole removes and returns the pole of type t closest to the unit circle.
func (pr *Pairing) takePole(t rootType) (p complex128) {
	p, pr.poles = takeAt(pr.poles, worstPole(pr.poles, t))
	return
}

// takeZero removes and returns the zero of type t nearest to target.
func (pr *Pairing) takeZero(target complex128, t rootType) (z complex128) {
	z, pr.zeros = takeAt(pr.zeros, nearestRoot(pr.zeros, target, t))
	return
}

// isReal reports whether r has no imaginary part, same as numpy.isreal.
// CollapseConjugates stores real roots with an imaginary part of exactly 0.
func isReal(r complex128) (ok bool) {
	ok = imag(r) == 0
	return
}

// countReal is how many of roots are real.
func countReal(roots []complex128) (n int) {
	for _, r := range roots {
		if isReal(r) {
			n++
		}
	}
	return
}

// worstPole is the index of the pole of type t closest to the unit circle,
// the first on ties, same as numpy.argmin.
func worstPole(poles []complex128, t rootType) (index int) {
	index = -1
	best := math.Inf(1)
	for i, p := range poles {
		if d := math.Abs(1 - cmplx.Abs(p)); t.matches(p) && d < best {
			best, index = d, i
		}
	}
	if index < 0 {
		panic(fmt.Sprintf("poly: no pole of type %d left to pair", t))
	}
	return
}

// nearestRoot is the index of the root of type t closest to target,
// the first on ties, same as numpy.argsort.
func nearestRoot(roots []complex128, target complex128, t rootType) (index int) {
	index = -1
	best := math.Inf(1)
	for i, r := range roots {
		if d := cmplx.Abs(r - target); t.matches(r) && d < best {
			best, index = d, i
		}
	}
	if index < 0 {
		panic(fmt.Sprintf("poly: no zero of type %d left to pair", t))
	}
	return
}

// takeAt removes and returns roots[i], keeping the rest in order, same as numpy.delete.
func takeAt(roots []complex128, i int) (r complex128, rest []complex128) {
	r = roots[i]
	rest = slices.Delete(roots, i, i+1)
	return
}

// padRoots returns a copy of roots extended to length n with roots at the origin.
func padRoots(roots []complex128, n int) (padded []complex128) {
	padded = make([]complex128, n)
	copy(padded, roots)
	return
}
