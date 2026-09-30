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
func butterworthPrototype(order int) (proto analogZPK) {
	proto.Poles = make([]complex128, order)
	for k := range order {
		m := float64(-order + 1 + 2*k) // -N+1, -N+3, ..., N-1
		theta := math.Pi * m / (2 * float64(order))
		proto.Poles[k] = -cmplx.Exp(complex(0, theta))
	}
	proto.Gain = 1
	return
}

// Butter designs a digital Butterworth filter for band and returns ZPK
// use its BA or SOS method for coefficients
//
// order must be >= 1
//
// Bandpass and Bandstop designs have order 2*order
// as each prototype pole becomes a pair
//
//	f, err := dsp.Butter(4, dsp.Lowpass(0.3))
//	f, a, err := dsp.Butter(2, dsp.Bandpass(0.2, 0.5))
//
// returns ErrInvalidOrder if order < 1, and ErrInvalidBand if band is zero-valued,
// has a frequency outside (0, 1), or has its low edge not below its high edge
//
// equivalent to scipy.signal.butter(order, wn, btype, output='zpk') and MATLAB's butter(order, wn, ftype)
func Butter(order int, band Band) (f ZPK, err error) {
	if order < 1 {
		return f, fmt.Errorf("%w: order must be >= 1, got %d", ErrInvalidOrder, order)
	}
	if err = band.validate(); err != nil {
		return f, err
	}
	f = band.transform(butterworthPrototype(order)).bilinear()
	return
}
