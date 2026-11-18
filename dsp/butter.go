package dsp

import (
	"math/cmplx"
)

// Butter designs a digital Butterworth filter for band:
// maximally flat in the passband, falling monotonically through the band edge,
// which is its -3 dB point.
// A Bandpass or Bandstop design has order 2·order, as each prototype pole becomes a pair.
// Use the ZPK's SOS method for filtering.
//
//	f, err := dsp.Butter(4, dsp.Lowpass(0.3))
//	sos, err := f.SOS()
//	y, err := dsp.SOSFilter(sos, x)
//
// Units:
//   - order is the prototype's order
//   - band's edges are fractions of Nyquist
//
// Returns:
//   - ErrInvalidOrder for an order below 1.
//   - ErrInvalidBand for the zero Band, an edge outside (0, 1), or a low edge not below the high.
//
// Equivalent to:
//   - scipy.signal.butter(order, wn, btype, output='zpk')
//   - MATLAB butter(order, wn, ftype)
func Butter(order int, band Band) (f ZPK, err error) {
	return designDigital(order, band, func() analogZPK { return butterworthPrototype(order) })
}

// butterworthPrototype returns the analog Butterworth lowpass prototype of the given order:
// -3 dB at 1 rad/s, no finite zeros, a gain of 1,
// and poles spaced evenly around the left half of the unit circle, same as scipy.signal.buttap.
func butterworthPrototype(order int) (proto analogZPK) {
	proto.Poles = make([]complex128, order)
	for k := range order {
		proto.Poles[k] = -cmplx.Exp(complex(0, prototypeAngle(order, k)))
	}
	proto.Gain = 1
	return
}
