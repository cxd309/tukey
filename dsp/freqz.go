package dsp

import (
	"fmt"
	"math"
	"math/cmplx"
)

// Freqz returns the frequency response of the filter b, a at n frequencies w,
// evenly spaced from 0 up to (not including) pi rad/sample, which is Nyquist
// h[k] is the complex gain at w[k], so |h[k]| is the magnitude and its angle the phase
//
// For frequencies normalised to Nyquist, as Band uses, divide w by pi
// for Hz at sample rate fs, multiply w by fs/(2*pi)
//
// returns ErrInvalidCoefficients as Filter does, and ErrInvalidLength if n < 0
//
// equivalent to scipy.signal.freqz(b, a, worN=n), and to MATLAB's [h, w] = freqz(b, a, n)
// (note MATLAB returns h first)
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

// SOSFreqz is Freqz for a filter given as second-order sections:
// the product of each section's response
//
// returns ErrInvalidSOS as SOSFilter does, and ErrInvalidLength if n < 0
//
// equivalent to scipy.signal.freqz_sos(sos, worN=n) (sosfreqz before SciPy 1.15)
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

// freqzPoints is n frequencies evenly spaced from 0 up to (not including) pi,
// as numpy.linspace(0, pi, n, endpoint=False)
func freqzPoints(n int) (w []float64) {
	w = make([]float64, n)
	step := math.Pi / float64(n)
	for k := range w {
		w[k] = float64(k) * step
	}
	return
}

// response evaluates B(z)/A(z) at z = e^(jw) for each frequency in w, with B and A
// as polynomials in z^-1 in the coefficients as given, as SciPy's freqz does for IIR filters
func response(b, a, w []float64) (h []complex128) {
	h = make([]complex128, len(w))
	for k, wk := range w {
		zInv := cmplx.Exp(complex(0, -wk))
		h[k] = horner(b, zInv) / horner(a, zInv)
	}
	return
}

// horner evaluates c[0] + c[1]x + c[2]x^2 + ... from the highest power down,
// as numpy.polynomial.polynomial.polyval does
func horner(c []float64, x complex128) (value complex128) {
	value = complex(c[len(c)-1], 0)
	for i := len(c) - 2; i >= 0; i-- {
		value = complex(c[i], 0) + value*x
	}
	return
}
