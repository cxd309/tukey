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
