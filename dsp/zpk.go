package dsp

import (
	"fmt"

	"github.com/cxd309/tukey/internal/poly"
)

// ZPK is a digital filter described by its zeros, poles and gain
// every IIR design function returns one or can be built by callers
// complex zeros and poles must come in conjugate pairs for the filter to be real
type ZPK struct {
	Zeros []complex128
	Poles []complex128
	Gain  float64
}

func (f ZPK) BA() (b, a []float64, err error) {
	bComplex := poly.FromRoots(f.Zeros)
	for i := range bComplex {
		bComplex[i] *= complex(f.Gain, 0)
	}
	if b, err = poly.RealCoeffs(bComplex); err != nil {
		return nil, nil, fmt.Errorf("%w: zeros: %v", ErrInvalidZPK, err)
	}
	if a, err = poly.RealCoeffs(poly.FromRoots(f.Poles)); err != nil {
		return nil, nil, fmt.Errorf("%w: poles: %v", ErrInvalidZPK, err)
	}
	return
}

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

	pr, sections, err := poly.NewPairing(f.Zeros, f.Poles)
	if err != nil {
		return nil, fmt.Errorf("%w: %v", ErrInvalidZPK, err)
	}

	sos = make(SOS, sections)
	// built last section first, so the poles closest to the unit circle end up last
	for si := sections - 1; si >= 0; si-- {
		zeros, poles := pr.Next()
		if sos[si], err = section(zeros, poles); err != nil {
			return nil, err
		}
	}

	for i := range 3 {
		sos[0][i] *= f.Gain
	}
	return
}

// section builds one second-order section from up to two zeros and two poles
// right-aligning each polynomial so a missing root leave a leading zero
// same as SciPy _single_zpksos
func section(zeros, poles []complex128) (s [6]float64, err error) {
	b, err := poly.RealCoeffs(poly.FromRoots(zeros))
	if err != nil {
		return s, fmt.Errorf("%w: %v", ErrInvalidZPK, err)
	}
	a, err := poly.RealCoeffs(poly.FromRoots(poles))
	if err != nil {
		return s, fmt.Errorf("%w: %v", ErrInvalidZPK, err)
	}
	copy(s[3-len(b):3], b)
	copy(s[6-len(a):6], a)
	return
}
