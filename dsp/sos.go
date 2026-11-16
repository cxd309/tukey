package dsp

import (
	"fmt"

	"github.com/cxd309/tukey/internal/poly"
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

// validate checks the sections as SciPy's sosfilt does:
// at least one section, and a0 == 1 in every one
func (s SOS) validate() (err error) {
	if len(s) == 0 {
		return fmt.Errorf("%w: no sections", ErrInvalidSOS)
	}
	for i, section := range s {
		if section[3] != 1 {
			return fmt.Errorf("%w: section %d has a0 = %v; every section's a0 must be 1", ErrInvalidSOS, i, section[3])
		}
	}
	return
}

// toCascade validates the sections as SciPy's sosfilt does and returns them as a cascade
func (s SOS) toCascade() (c cascade, err error) {
	if err = s.validate(); err != nil {
		return nil, err
	}
	c = make(cascade, len(s))
	for i, section := range s {
		if c[i], err = newDigitalFilter(section[:3], section[3:]); err != nil {
			return nil, err
		}
	}
	return
}

// ntaps is SciPy's tap count for sosfiltfilt's default padding:
// 2*sections + 1 minus 1 for each first-order section (b2 and a2 both zero)
func (s SOS) ntaps() (n int) {
	zeroB2, zeroA2 := 0, 0
	for _, section := range s {
		if section[2] == 0 {
			zeroB2++
		}
		if section[5] == 0 {
			zeroA2++
		}
	}
	n = 2*len(s) + 1 - min(zeroB2, zeroA2)
	return
}

// size is the filter's size for its padding length:
// scipy ntaps, and matlab order, summed over the sections
func (s SOS) size() (size filterSize) {
	size.ntaps = s.ntaps()
	for _, section := range s {
		size.order += effectiveOrder(section[:3], section[3:])
	}
	return
}
