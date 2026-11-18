package dsp

import (
	"fmt"
	"math"
	"math/cmplx"

	"github.com/cxd309/tukey/internal/poly"
)

// ZPK is a digital filter described by its zeros, poles and gain:
//
//	H(z) = Gain·Π (z - zero)/Π (z - pole)
//
// Every IIR design returns one, and one can be built by hand.
// Complex zeros and poles must come in conjugate pairs, for the filter to be real.
// Use its SOS method for filtering, and its BA method for low orders or for interchange.
type ZPK struct {
	Zeros []complex128 // in the z-plane
	Poles []complex128 // in the z-plane, inside the unit circle for a stable filter
	Gain  float64      // dimensionless
}

// BA returns the filter's transfer-function coefficients, b over a, from its zeros and poles,
// with a[0] = 1.
// As the order rises or the band narrows, b and a lose digits, and SOS keeps them.
//
//	b, a, err := f.BA()
//	y, err := dsp.Filter(b, a, x)
//
// Units: b and a are dimensionless
//
// Returns: ErrInvalidZPK if complex zeros or poles aren't in conjugate pairs, or the gain, a zero or a pole isn't finite.
//
// Equivalent to:
//   - scipy.signal.zpk2tf(z, p, k)
//   - MATLAB zp2tf(z, p, k)
func (f ZPK) BA() (b, a []float64, err error) {
	if err = f.validate(); err != nil {
		return nil, nil, err
	}
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

// SOS returns the filter as second-order sections, for filtering at any order.
// Each pole is paired with its nearest zero, starting from the poles closest to the unit circle,
// which go in the last sections.
// That keeps each section's peak gain small,
// so SOS stays accurate at orders where BA's coefficients don't.
// The algorithm follows SciPy's, as THIRD_PARTY_NOTICES records.
//
//	f, err := dsp.Butter(8, dsp.Lowpass(0.1))
//	sos, err := f.SOS()
//	y, err := dsp.SOSFilter(sos, x)
//
// Units: sos's coefficients are dimensionless
//
// Returns: ErrInvalidZPK if complex zeros or poles aren't in conjugate pairs, or the gain, a zero or a pole isn't finite.
//
// Equivalent to:
//   - scipy.signal.zpk2sos(z, p, k), with its default pairing='nearest'
//   - MATLAB zp2sos(z, p, k)
func (f ZPK) SOS() (sos SOS, err error) {
	if err = f.validate(); err != nil {
		return nil, err
	}
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

// section builds one second-order section from up to two zeros and two poles,
// right-aligning each polynomial so a missing root leaves a leading 0, same as SciPy's _single_zpksos.
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

// validate returns ErrInvalidZPK unless the gain, zeros and poles are all finite.
func (f ZPK) validate() (err error) {
	if math.IsNaN(f.Gain) || math.IsInf(f.Gain, 0) {
		return fmt.Errorf("%w: gain %g, must be finite", ErrInvalidZPK, f.Gain)
	}
	for _, roots := range [][]complex128{f.Zeros, f.Poles} {
		for _, r := range roots {
			if cmplx.IsNaN(r) || cmplx.IsInf(r) {
				return fmt.Errorf("%w: root %v, must be finite", ErrInvalidZPK, r)
			}
		}
	}
	return
}
