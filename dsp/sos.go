package dsp

import (
	"fmt"
)

// SOS is a filter as a cascade of second-order sections, applied in order.
// Each row is one section's coefficients, b0, b1, b2, a0, a1, a2, the same layout as SciPy's sos arrays,
// with a0 = 1.
// Its coefficients are dimensionless.
// ZPK's SOS method builds one from a design.
type SOS [][6]float64

// validate returns ErrInvalidSOS unless there's at least one section,
// and every section has a0 = 1 and is finite.
// SciPy's sosfilt checks the first two.
func (s SOS) validate() (err error) {
	if len(s) == 0 {
		return fmt.Errorf("%w: no sections", ErrInvalidSOS)
	}
	for i, section := range s {
		if section[3] != 1 {
			return fmt.Errorf("%w: section %d has a0 = %v, and every section's a0 must be 1", ErrInvalidSOS, i, section[3])
		}
		if !allFinite(section[:]) {
			return fmt.Errorf("%w: section %d must be finite", ErrInvalidSOS, i)
		}
	}
	return
}

// toCascade validates the sections and returns them as a cascade.
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
// 2·sections + 1, less 1 for each first-order section, whose b2 and a2 are both 0.
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
// SciPy's ntaps, and MATLAB's order summed over the sections.
func (s SOS) size() (size filterSize) {
	size.ntaps = s.ntaps()
	for _, section := range s {
		size.order += effectiveOrder(section[:3], section[3:])
	}
	return
}
