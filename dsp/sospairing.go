package dsp

import "math/cmplx"

// pairing holds the zeros and poles not yet placed in a section, with each
// conjugate pair stored once (see collapseConjugates); follows the loop in
// SciPy's zpk2sos
type pairing struct {
	zeros, poles []complex128
}

// newPairing balances f's zeros and poles and collapses conjugate pairs,
// returning the pool and how many sections it makes
func newPairing(f ZPK) (pr pairing, sections int, err error) {
	// balance the counts with roots at the origin; an odd count gets one more of
	// each, so every section has exactly two zeros and two poles
	n := max(len(f.Zeros), len(f.Poles))
	zeros, poles := padRoots(f.Zeros, n), padRoots(f.Poles, n)
	sections = (n + 1) / 2
	if n%2 == 1 {
		zeros, poles = append(zeros, 0), append(poles, 0)
	}

	if pr.zeros, err = collapseConjugates(zeros); err != nil {
		return pr, 0, err
	}
	pr.poles, err = collapseConjugates(poles)
	return
}

// next removes and returns the zeros and poles for the next section, starting
// from the pole closest to the unit circle
func (pr *pairing) next() (zeros, poles []complex128) {
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

// lastRealPole pairs the only remaining real pole with the nearest real zero,
// plus a pole and zero at the origin
func (pr *pairing) lastRealPole(p1 complex128) (zeros, poles []complex128) {
	zeros = []complex128{pr.takeZero(p1, rootReal), 0}
	poles = []complex128{p1, 0}
	return
}

// mustKeepRealPairTogether reports whether exactly one real pole and one real zero
// remain for a later section, in which case complex pole p1 must take a complex zero
func (pr *pairing) mustKeepRealPairTogether(p1 complex128) (must bool) {
	must = !isReal(p1) && len(pr.poles)+1 == len(pr.zeros) &&
		countReal(pr.poles) == 1 && countReal(pr.zeros) == 1
	return
}

// complexPair pairs complex pole p1 with the nearest complex zero, and both with their conjugates
func (pr *pairing) complexPair(p1 complex128) (zeros, poles []complex128) {
	z1 := pr.takeZero(p1, rootComplex)
	zeros = []complex128{z1, cmplx.Conj(z1)}
	poles = []complex128{p1, cmplx.Conj(p1)}
	return
}

// secondOrder completes a section around p1: its conjugate (or, if real, the next
// real pole closest to the unit circle), then the nearest zero and whatever
// completes it: its conjugate if complex, or the nearest other real zero
func (pr *pairing) secondOrder(p1 complex128) (zeros, poles []complex128) {
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

// takePole removes and returns the pole of type t closest to the unit circle
func (pr *pairing) takePole(t rootType) (p complex128) {
	p, pr.poles = takeAt(pr.poles, worstPole(pr.poles, t))
	return
}

// takeZero removes and returns the zero of type t nearest to target
func (pr *pairing) takeZero(target complex128, t rootType) (z complex128) {
	z, pr.zeros = takeAt(pr.zeros, nearestRoot(pr.zeros, target, t))
	return
}
