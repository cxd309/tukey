package dsp

import (
	"fmt"
	"math"
	"math/cmplx"
)

// Freqz returns the frequency response of the filter b/a at n frequencies w,
// evenly spaced from 0 up to, but not including, Nyquist.
// h[k] is the complex gain at w[k], so |h[k]| is the magnitude, and its angle the phase.
// n = 0 gives empty slices, same as SciPy, where MATLAB's freqz errors.
//
//	w, h, err := dsp.Freqz(b, a, 512)
//	f, err := dsp.RadPerSampleToHz(w, fs)
//
// Units:
//   - b and a are dimensionless
//   - w is in rad/sample, where π is Nyquist, and RadPerSampleToHz or RadPerSampleToNyquist converts it
//   - h is dimensionless, a gain
//
// Returns:
//   - ErrInvalidCoefficients for an empty b or a, an a[0] of 0, or a NaN or infinite coefficient.
//   - ErrInvalidLength for a negative n.
//
// Equivalent to:
//   - scipy.signal.freqz(b, a, worN=n)
//   - MATLAB [h, w] = freqz(b, a, n), which returns h first
func Freqz(b, a []float64, n int) (w []float64, h []complex128, err error) {
	if err = validateCoefficients(b, a); err != nil {
		return nil, nil, err
	}
	if n < 0 {
		return nil, nil, fmt.Errorf("%w: n must be >= 0, got %d", ErrInvalidLength, n)
	}
	w = freqzPoints(n)
	h = response(b, a, w)
	return
}

// SOSFreqz returns the frequency response of the filter sos at n frequencies w,
// evenly spaced from 0 up to, but not including, Nyquist:
// the product of each section's response, which stays accurate where a high order's (b, a) doesn't.
// n = 0 gives empty slices.
//
//	w, h, err := dsp.SOSFreqz(sos, 512)
//
// Units:
//   - sos is dimensionless
//   - w is in rad/sample, where π is Nyquist, and RadPerSampleToHz or RadPerSampleToNyquist converts it
//   - h is dimensionless, a gain
//
// Returns:
//   - ErrInvalidSOS for no sections, a section whose a0 isn't 1, or a NaN or infinite coefficient.
//   - ErrInvalidLength for a negative n.
//
// Equivalent to:
//   - scipy.signal.freqz_sos(sos, worN=n), called sosfreqz before SciPy 1.15
//   - MATLAB [h, w] = freqz(sos, n), which returns h first
func SOSFreqz(sos SOS, n int) (w []float64, h []complex128, err error) {
	if err = sos.validate(); err != nil {
		return nil, nil, err
	}
	if n < 0 {
		return nil, nil, fmt.Errorf("%w: n must be >= 0, got %d", ErrInvalidLength, n)
	}
	w = freqzPoints(n)
	h = make([]complex128, n)
	for k := range h {
		h[k] = 1
	}
	for _, section := range sos {
		sectionH := response(section[:3], section[3:], w)
		for k := range h {
			h[k] *= sectionH[k]
		}
	}
	return
}

// freqzPoints is n frequencies evenly spaced from 0 up to, but not including, π,
// same as numpy.linspace(0, pi, n, endpoint=False).
func freqzPoints(n int) (w []float64) {
	w = make([]float64, n)
	step := math.Pi / float64(n)
	for k := range w {
		w[k] = float64(k) * step
	}
	return
}

// response evaluates B(z)/A(z) at z = e^(jw) for each frequency in w,
// with B and A as polynomials in z^-1, same as SciPy's freqz for an IIR filter.
func response(b, a, w []float64) (h []complex128) {
	h = make([]complex128, len(w))
	for k, wk := range w {
		zInv := cmplx.Exp(complex(0, -wk))
		h[k] = horner(b, zInv) / horner(a, zInv)
	}
	return
}

// horner evaluates c[0] + c[1]·x + c[2]·x² + ... from the highest power down,
// same as numpy.polynomial.polynomial.polyval.
func horner(c []float64, x complex128) (value complex128) {
	value = complex(c[len(c)-1], 0)
	for i := len(c) - 2; i >= 0; i-- {
		value = complex(c[i], 0) + value*x
	}
	return
}
