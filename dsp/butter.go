package dsp

import (
	"fmt"
	"math"
	"math/cmplx"
)

// butterworthPrototype returns the analog Butterworth lowpass prototype of the given order:
// -3dB at 1 rad/s, no finite zeros, unity gain,
// and poles spaced evenly around the left half of the unit circle
// equivalent to scipy.signal.buttap
func butterworthPrototype(order int) (proto zpk) {
	proto.poles = make([]complex128, order)
	for k := range order {
		m := float64(-order + 1 + 2*k) // -N+1, -N+3, ..., N-1
		theta := math.Pi * m / (2 * float64(order))
		proto.poles[k] = -cmplx.Exp(complex(0, theta))
	}
	proto.gain = 1
	return
}

// Butter designs a digital Butterworth filter for band and returns its
// transfer-function coefficients b (numerator) and a (denominator)
//
// order must be >= 1
//
// Bandpass and Bandstop designs have order 2*order
// as each prototype pole becomes a pair
//
//	b, a, err := dsp.Butter(4, dsp.Lowpass(0.3))
//	b, a, err := dsp.Butter(2, dsp.Bandpass(0.2, 0.5))
//
// returns ErrInvalidOrder if order < 1, and ErrInvalidBand if band is zero-valued,
// has a frequency outside (0, 1), or has its low edge not below its high edge
//
// equivalent to scipy.signal.butter(order, wn, btype) and MATLAB's butter(order, wn, ftype)
func Butter(order int, band Band) (b, a []float64, err error) {
	if order < 1 {
		return nil, nil, fmt.Errorf("%w: order must be >= 1, got %d", ErrInvalidOrder, order)
	}
	if err = band.validate(); err != nil {
		return nil, nil, err
	}

	analog := band.transform(butterworthPrototype(order))
	b, a = analog.bilinear().transferFunction()
	return
}
