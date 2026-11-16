package dsp

import (
	"fmt"
	"math"

	"github.com/cxd309/tukey/internal/floats"
	"github.com/cxd309/tukey/internal/special"
)

// FIRWin designs a linear-phase FIR filter for band by the window method and
// returns its taps: the ideal response's impulse response, truncated to the
// window's length and multiplied by it.
// The window sets the transition width and stopband attenuation,
// and its length, the number of taps, the sharpness:
//
//	w, err := window.Hamming(101)
//	h, err := dsp.FIRWin(w, dsp.Lowpass(0.2))
//	y, err := dsp.Convolve(x, h, dsp.ConvolveSame())
//
// The taps are scaled so the response is exactly 1 at DC if the band passes
// it, at Nyquist for a highpass, and otherwise at the passband's centre.
// FIRNoScale turns that off.
// Any window is accepted, but only a symmetric one gives exactly linear phase:
// build it with window's default symmetric form, not Periodic.
// Its delay is (len(window) − 1)/2 samples.
//
// Returns:
//   - ErrInvalidLength if window is empty, or has an even length for a
//     Highpass or Bandstop, whose response at Nyquist an even-length filter
//     can't pass
//   - ErrInvalidBand for a band Butter rejects
//
// Equivalent to:
//   - scipy.signal.firwin(len(window), wn, window, pass_zero=btype), whose
//     window is given by name
//   - matlab fir1(len(window) − 1, wn, ftype, window)
func FIRWin(window []float64, band Band, opts ...FIROption) (h []float64, err error) {
	var cfg firConfig
	for _, opt := range opts {
		opt(&cfg)
	}
	if err = band.validate(); err != nil {
		return nil, err
	}
	n := len(window)
	passesNyquist := band.btype == highpassType || band.btype == bandstopType
	switch {
	case n == 0:
		return nil, fmt.Errorf("%w: FIRWin needs a window of at least 1 sample", ErrInvalidLength)
	case passesNyquist && n%2 == 0:
		return nil, fmt.Errorf("%w: FIRWin: %v needs an odd number of taps, got %d", ErrInvalidLength, band, n)
	}

	bands, scaleAt := band.passbands()
	h = make([]float64, n)
	centre := float64(n-1) / 2
	for k, w := range window {
		m := float64(k) - centre // exact: a whole or half number
		var ideal float64
		for _, p := range bands {
			ideal += p.response(m)
		}
		h[k] = w * ideal
	}
	if cfg.noScale {
		return
	}

	// the zero-phase amplitude at scaleAt, Σ h·cos(π·m·scaleAt)
	var amplitude floats.Sum
	for k, tap := range h {
		_, cos := special.SinCosPiProduct(scaleAt, float64(k)-centre)
		amplitude.Add(float64(tap * cos))
	}
	scale := amplitude.Total()
	for k := range h {
		h[k] /= scale
	}
	return
}

// FIROption configures FIRWin
type FIROption func(*firConfig)

// FIRNoScale leaves FIRWin's taps unscaled: the windowed ideal response as it is
func FIRNoScale() FIROption {
	return func(c *firConfig) { c.noScale = true }
}

// firConfig is FIRWin's settings, built from FIROptions
type firConfig struct {
	noScale bool
}

// passbands is the bands band passes, and where a scaled design has gain 1:
// DC if the first reaches it, Nyquist if it reaches that, otherwise its centre
func (band Band) passbands() (bands []passband, scaleAt float64) {
	edges := band.edges
	switch band.btype {
	case lowpassType:
		bands = []passband{{0, edges[0]}}
	case highpassType:
		bands = []passband{{edges[0], 1}}
	case bandpassType:
		bands = []passband{{edges[0], edges[1]}}
	default:
		bands = []passband{{0, edges[0]}, {edges[1], 1}}
	}
	switch first := bands[0]; {
	case first.left == 0:
		scaleAt = 0
	case first.right == 1:
		scaleAt = 1
	default:
		scaleAt = (first.left + first.right) / 2
	}
	return
}

// passband is one band an ideal FIR response passes, from left to right,
// normalised to Nyquist
type passband struct {
	left, right float64
}

// response is the ideal response's impulse response at m samples from the centre:
// for each passband, (sin(π·right·m) − sin(π·left·m))/(π·m), or right − left at m = 0.
// The difference of sines is formed as the product 2·cos(π·centre·m)·sin(π·halfWidth·m),
// which doesn't cancel for a narrow band.
func (p passband) response(m float64) (h float64) {
	if m == 0 {
		return p.right - p.left
	}
	_, cos := special.SinCosPiProduct((p.left+p.right)/2, m)
	sin, _ := special.SinCosPiProduct((p.right-p.left)/2, m)
	return 2 * cos * sin / (math.Pi * m)
}
