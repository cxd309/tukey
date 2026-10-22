package dsp

import (
	"fmt"
	"math"
	"math/cmplx"

	"github.com/cxd309/tukey/internal/floats"
)

// EnvelopeMethod is how Envelope measures a signal's envelope:
// EnvelopeAnalytic (the zero EnvelopeMethod) or EnvelopeRMS.
type EnvelopeMethod struct {
	rms    bool
	window int
}

// EnvelopeAnalytic is the magnitude of the analytic signal (see Hilbert).
// It follows each scillation's peaks closely, so suits a signal with one main frequency
func EnvelopeAnalytic() (method EnvelopeMethod) {
	return
}

// EnvelopeRMS is the moving RMS over a window of n samples centred on each sample.
// For an even n the extra sample is after, so n = 10 covers k-4 to k+5.
// The window shrinks where it would pass either end of the signal.
func EnvelopeRMS(n int) (method EnvelopeMethod) {
	method = EnvelopeMethod{rms: true, window: n}
	return
}

// Envelope returns the upper and lower enveopes of x:
// its mean plus and minus the envelope of x less its mean, measured according method.
//
// Returns:
//   - ErrInvalidLength for an empty x
//   - ErrInvalidParameter for an EnvelopeRMS window below1
//
// Equivalent to MATLAB [upper, lower] = envelope(x) and envelope(x, n, "rms").
// scipy.signal.envelope drops the Nyquist bin from the analytic envelope, and has no RMS method.
func Envelope(x []float64, method EnvelopeMethod) (upper, lower []float64, err error) {
	if len(x) == 0 {
		return nil, nil, fmt.Errorf("%w: Envelope of an empty signal", ErrInvalidLength)
	}
	if method.rms && method.window < 1 {
		return nil, nil, fmt.Errorf("%w: EnvelopeRMS window %d, must be at least 1", ErrInvalidParameter, method.window)
	}
	m := mean(x)
	centred := make([]float64, len(x))
	for i, v := range x {
		centred[i] = v - m
	}
	var env []float64
	if method.rms {
		env = movingRMS(centred, method.window)
	} else if env, err = analyticEnvelope(centred); err != nil {
		return nil, nil, err
	}
	upper = make([]float64, len(x))
	lower = make([]float64, len(x))
	for i, e := range env {
		upper[i], lower[i] = m+e, m-e
	}
	return
}

// analyticEnvelope is the magnitude of x's analytic signal.
func analyticEnvelope(x []float64) (env []float64, err error) {
	analytic, err := Hilbert(x)
	if err != nil {
		return nil, err
	}
	env = make([]float64, len(x))
	for i, a := range analytic {
		env[i] = cmplx.Abs(a)
	}
	return
}

// movingRMS is the RMS of x over a window of n samples around each one,
// (n-1)/2 before and n/2 after, shrinking where it passes either end of x.
// Each window's squares are summed with compensation, so one large square
// doesn't loose the low bits in small squares
func movingRMS(x []float64, n int) (rms []float64) {
	before, after := (n-1)/2, n/2
	rms = make([]float64, len(x))
	for k := range x {
		lo, hi := max(0, k-before), min(len(x), k+after+1)
		var sum floats.Sum
		for _, v := range x[lo:hi] {
			sum.Add(v * v)
		}
		rms[k] = math.Sqrt(sum.Total() / float64(hi-lo))
	}
	return
}
