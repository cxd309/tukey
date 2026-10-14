package fft

import (
	"fmt"
	"math"
)

// Normalisation chooses where a tranform pair's 1/n scaling goes
type Normalisation int

const (
	// NormBackward leaves the forward tranform unscaled and scales the inverse by 1/n
	// the default (same as SciPy and MATLAB)
	NormBackward Normalisation = iota
	// NormOrtho scales both by 1/sqrt(n)
	// unitary tranforms and same signal energy in both domains
	NormOrtho
	// NormForward scales the forward tranform by 1/n
	// so DC bin is signal's mean and leave inverse unchanged
	NormForward
)

// Option configures a tranform (see Norm and N)
type Option func(*config)

// config is a tranform's settings
// built from Options
type config struct {
	norm Normalisation
	n    int  // tranform length
	nSet bool // if N was given, otherwise n is len(input)
}

// Norm sets the normalisation, default is NormBackward
func Norm(norm Normalisation) (option Option) {
	option = func(c *config) { c.norm = norm }
	return
}

// N sets the transform length, input is truncated or zero-padded to n first
// not accepted by InverseReal as length is its own n argument
// equivalent to SciPy n argument and MATLAB fft(x, n)
func N(n int) (option Option) {
	option = func(c *config) { c.n, c.nSet = n, true }
	return
}

func newConfig(opts []Option, inputLength int) (cfg config, err error) {
	for _, opt := range opts {
		opt(&cfg)
	}
	if !(cfg.norm >= NormBackward && cfg.norm <= NormForward) {
		return cfg, fmt.Errorf("%w: unknown normalisation %d", ErrInvalidOption, int(cfg.norm))
	}
	if !cfg.nSet {
		cfg.n = inputLength
	}
	if cfg.n < 1 {
		return cfg, fmt.Errorf("%w: tranform length must be >= 1, got %d", ErrInvalidLength, cfg.n)
	}
	return
}

// scale is the factor a tranform of length n multiplies its output by in given direction under this norm
func (norm Normalisation) scale(n int, forward bool) (s float64) {
	switch {
	case norm == NormOrtho:
		s = 1 / math.Sqrt(float64(n))
	case (norm == NormBackward) != forward: // backward-scaled inverse, or forward-scaled forward
		s = 1 / float64(n)
	default:
		s = 1
	}
	return
}

// scaled multiplies each value in y by s in place
// exactly nothing changes when s is 1
func scaled(y []complex128, s float64) {
	if s == 1 {
		return
	}
	for i := range y {
		y[i] *= complex(s, 0) // multiplied, not divded, same as pocketfft
	}
}

// resized returns a copy of x truncated or zero-padded to length n
func resized[T float64 | complex128](x []T, n int) (y []T) {
	y = make([]T, n)
	copy(y, x)
	return
}
