package dsp

import (
	"fmt"
	"math"
)

// IIRNotch designs a second-order notch filter,
// which removes the frequency w0 and passes the rest at 0 dB.
// The quality factor q sets the -3 dB bandwidth, w0/q, so the higher q, the narrower the notch.
// Its zeros are on the unit circle at w0, and its poles just inside them.
// For mains hum at 50 Hz, sampled at 1 kHz:
//
//	w0, err := dsp.HzToNyquist(50.0, 1000)
//	f, err := dsp.IIRNotch(w0, 30)
//
// Units:
//   - w0 is a fraction of Nyquist
//   - q is dimensionless
//
// Returns: ErrInvalidParameter for a w0 outside (0, 1), a q that isn't positive and finite, or a bandwidth w0/q that isn't below Nyquist.
//
// Equivalent to:
//   - scipy.signal.iirnotch(w0, q), which returns (b, a) where IIRNotch returns zeros, poles and gain
//   - MATLAB iirnotch(w0, w0/q), in the DSP System Toolbox
func IIRNotch(w0, q float64) (f ZPK, err error) {
	d, err := newNotchDesign("IIRNotch", w0, q)
	if err != nil {
		return f, err
	}
	f = ZPK{
		Zeros: []complex128{complex(d.cos, d.sin), complex(d.cos, -d.sin)},
		Poles: d.poles(),
		Gain:  1 / (1 + d.beta),
	}
	return
}

// IIRPeak designs a second-order peak filter,
// which passes the frequency w0 at 0 dB and rejects DC and Nyquist.
// The quality factor q sets the -3 dB bandwidth, w0/q, so the higher q, the narrower the peak.
//
//	f, err := dsp.IIRPeak(0.3, 10)
//
// Units:
//   - w0 is a fraction of Nyquist
//   - q is dimensionless
//
// Returns: ErrInvalidParameter for a w0 outside (0, 1), a q that isn't positive and finite, or a bandwidth w0/q that isn't below Nyquist.
//
// Equivalent to:
//   - scipy.signal.iirpeak(w0, q), which returns (b, a), and forms the gain as 1 - 1/(1 + β),
//     losing digits for a narrow bandwidth, 6e-13 at w0 = 0.002, q = 30.
//     IIRPeak's is β/(1 + β), exact to rounding.
//   - MATLAB iirpeak(w0, w0/q), in the DSP System Toolbox
func IIRPeak(w0, q float64) (f ZPK, err error) {
	d, err := newNotchDesign("IIRPeak", w0, q)
	if err != nil {
		return f, err
	}
	f = ZPK{
		Zeros: []complex128{1, -1},
		Poles: d.poles(),
		Gain:  d.beta / (1 + d.beta),
	}
	return
}

// newNotchDesign checks w0 and q, returning ErrInvalidParameter for the function called name,
// and returns their design.
func newNotchDesign(name string, w0, q float64) (d notchDesign, err error) {
	switch {
	case !(w0 > 0 && w0 < 1):
		err = fmt.Errorf("%w: %s: w0 must be in (0, 1), normalised to Nyquist, got %g", ErrInvalidParameter, name, w0)
	case !(q > 0) || math.IsInf(q, 1):
		err = fmt.Errorf("%w: %s: q must be positive and finite, got %g", ErrInvalidParameter, name, q)
	case !(w0/q < 1):
		err = fmt.Errorf("%w: %s: the bandwidth w0/q must be below Nyquist, got %g/%g", ErrInvalidParameter, name, w0, q)
	}
	if err != nil {
		return d, err
	}
	d.sin, d.cos = math.Sincos(math.Pi * w0)
	d.beta = math.Tan(math.Pi * w0 / (2 * q))
	return
}

// notchDesign is what IIRNotch and IIRPeak share:
// the centre and bandwidth, and the poles they put there.
type notchDesign struct {
	cos, sin float64 // of π·w0, the centre on the unit circle
	beta     float64 // tan(π·w0/(2q)), from the -3 dB bandwidth w0/q
}

// poles are the roots of z² - 2g·cos·z + (1 - β)/(1 + β), for g = 1/(1 + β).
// Their discriminant, g²·cos² - (1 - β)/(1 + β), is (β² - sin²)/(1 + β)²,
// so a narrow notch's poles, hard against the unit circle,
// come from (sin - β)(sin + β) rather than a difference of nearly equal squares.
// A bandwidth wide enough that β > sin gives two real poles,
// the smaller from their product, (1 - β)/(1 + β), rather than by cancellation.
func (d notchDesign) poles() (poles []complex128) {
	scale := 1 + d.beta
	if d.beta < d.sin {
		re := d.cos / scale
		im := math.Sqrt((d.sin-d.beta)*(d.sin+d.beta)) / scale
		return []complex128{complex(re, im), complex(re, -im)}
	}
	larger := (d.cos + math.Copysign(math.Sqrt((d.beta-d.sin)*(d.beta+d.sin)), d.cos)) / scale
	smaller := (1 - d.beta) / scale / larger
	return []complex128{complex(larger, 0), complex(smaller, 0)}
}
