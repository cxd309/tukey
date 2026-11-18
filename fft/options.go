package fft

import (
	"fmt"
	"math"
)

// Option configures a transform: NormBackward, NormOrtho, NormForward or N.
type Option func(*config)

// NormBackward leaves the forward transform unscaled and scales the inverse by 1/n.
// It's the default.
//
// Equivalent to:
//   - scipy.fft's norm="backward", its default
//   - MATLAB fft and ifft, which have no other
func NormBackward() (option Option) {
	option = func(c *config) { c.norm = normBackward }
	return
}

// NormOrtho scales both transforms by 1/√n,
// so they're unitary and the signal's energy is the same in both domains.
//
// Equivalent to:
//   - scipy.fft's norm="ortho"
//   - MATLAB fft(x)/sqrt(n)
func NormOrtho() (option Option) {
	option = func(c *config) { c.norm = normOrtho }
	return
}

// NormForward scales the forward transform by 1/n and leaves the inverse unscaled,
// so the DC bin is the signal's mean.
//
// Equivalent to:
//   - scipy.fft's norm="forward"
//   - MATLAB fft(x)/n
func NormForward() (option Option) {
	option = func(c *config) { c.norm = normForward }
	return
}

// N sets the transform length, truncating or zero-padding the input to n samples first.
// InverseReal doesn't take it, as its length is its n argument.
//
// Units: n is a number of samples
//
// Equivalent to:
//   - scipy.fft.fft(x, n)
//   - MATLAB fft(x, n)
func N(n int) (option Option) {
	option = func(c *config) { c.n, c.nSet = n, true }
	return
}

// normalisation is where a transform pair's 1/n scaling goes.
type normalisation int

const (
	normBackward normalisation = iota
	normOrtho
	normForward
)

// config is a transform's settings, built from Options.
type config struct {
	norm normalisation
	n    int  // transform length
	nSet bool // whether N was given, otherwise n is the input's length
}

// newConfig builds a transform's settings from opts,
// returning ErrInvalidLength for an N below 1, and ErrSignalTooShort for an empty input without N.
func newConfig(opts []Option, inputLength int) (cfg config, err error) {
	for _, opt := range opts {
		opt(&cfg)
	}
	switch {
	case cfg.nSet && cfg.n < 1:
		return cfg, fmt.Errorf("%w: N %d, must be at least 1", ErrInvalidLength, cfg.n)
	case !cfg.nSet && inputLength < 1:
		return cfg, fmt.Errorf("%w: cannot transform an empty signal", ErrSignalTooShort)
	case !cfg.nSet:
		cfg.n = inputLength
	}
	return
}

// scale is the factor a transform of length n multiplies its output by, in the given direction.
func (norm normalisation) scale(n int, forward bool) (s float64) {
	switch {
	case norm == normOrtho:
		s = 1 / math.Sqrt(float64(n))
	case (norm == normBackward) != forward: // backward-scaled inverse, or forward-scaled forward
		s = 1 / float64(n)
	default:
		s = 1
	}
	return
}

// scaled multiplies each value in y by s, in place.
// Nothing changes when s is 1.
func scaled(y []complex128, s float64) {
	if s == 1 {
		return
	}
	for i := range y {
		y[i] *= complex(s, 0) // multiplied, not divided, same as pocketfft
	}
}

// resized returns a copy of x, truncated or zero-padded to length n.
func resized[T float64 | complex128](x []T, n int) (y []T) {
	y = make([]T, n)
	copy(y, x)
	return
}
