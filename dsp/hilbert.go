package dsp

import (
	"fmt"

	"github.com/cxd309/tukey/fft"
)

// Hilbert returns the analytic signal of x:
// x + (i times x's hilbert transform)
// which is x frequency shifted by -90 degrees.
// Its magnitude is x's envelope and its angle x's instantaneous phase.
//
// It works in the frequency domain, doubling the positive frequencies and
// zeroing the negativeones, so it treats x as one period of a periodic signal.
//
// Returns ErrInvalidLength for an empty x.
//
// Equivalent to scipy.signal.hilbert(x), MATLAB hilbert(x)
func Hilbert(x []float64) (analytic []complex128, err error) {
	n := len(x)
	if n == 0 {
		return nil, fmt.Errorf("%w: Hilbert of an empty signal", ErrInvalidLength)
	}
	half, err := fft.ForwardReal(x) // bins 0 to n/2
	if err != nil {
		return nil, err
	}
	spectrum := make([]complex128, n) // the negative frequencies stay 0
	spectrum[0] = half[0]
	for k := 1; k < (n+1)/2; k++ {
		spectrum[k] = 2 * half[k]
	}
	if n%2 == 0 {
		spectrum[n/2] = half[n/2]
	} // Nyquist is its own negative frequency, so it's kept not doubled}
	analytic, err = fft.Inverse(spectrum)
	return
}
