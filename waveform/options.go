package waveform

import (
	"fmt"
	"strings"
)

// Option configures a waveform.
// One Option type serves every waveform, and each option means the same wherever it applies.
// A waveform given an option that doesn't apply to it returns ErrInvalidOption.
type Option func(*config)

// Linear sweeps Chirp's frequency linearly from f0 at t = 0 to f1 at t = t1:
// f(t) = f0 + (f1 - f0)·t/t1.
// It's the default.
//
// Equivalent to:
//   - scipy.signal.chirp's method='linear', its default
//   - MATLAB chirp's "linear", its default
func Linear() (option Option) {
	option = func(c *config) {
		c.sweep = sweepLinear
		c.given[optionLinear] = true
	}
	return
}

// Quadratic sweeps Chirp's frequency along a parabola from f0 at t = 0 to f1 at t = t1,
// with its vertex at t = 0, f(t) = f0 + (f1 - f0)·(t/t1)²,
// or with VertexAtT1 at t = t1, f(t) = f1 - (f1 - f0)·(1 - t/t1)².
//
// Equivalent to:
//   - scipy.signal.chirp's method='quadratic'
//   - MATLAB chirp's "quadratic", with its default shape, for the vertex at t = 0
func Quadratic() (option Option) {
	option = func(c *config) {
		c.sweep = sweepQuadratic
		c.given[optionQuadratic] = true
	}
	return
}

// Logarithmic sweeps Chirp's frequency geometrically from f0 at t = 0 to f1 at t = t1:
// f(t) = f0·(f1/f0)^(t/t1).
// f0 and f1 must be non-zero and of the same sign.
//
// Equivalent to:
//   - scipy.signal.chirp's method='logarithmic'
//   - MATLAB chirp's "logarithmic", except at f0 = f1, where MATLAB divides by 0
func Logarithmic() (option Option) {
	option = func(c *config) {
		c.sweep = sweepLogarithmic
		c.given[optionLogarithmic] = true
	}
	return
}

// Hyperbolic sweeps Chirp's frequency from f0 at t = 0 to f1 at t = t1 with a constant change in period:
// f(t) = f0·f1·t1/((f0 - f1)·t + f1·t1).
// f0 and f1 must be non-zero.
// The frequency is infinite at t = f1·t1/(f1 - f0), which for an upward sweep is after t1.
//
// Equivalent to:
//   - scipy.signal.chirp's method='hyperbolic'
//   - MATLAB has none
func Hyperbolic() (option Option) {
	option = func(c *config) {
		c.sweep = sweepHyperbolic
		c.given[optionHyperbolic] = true
	}
	return
}

// VertexAtT1 puts a Quadratic sweep's vertex at t = t1, where the frequency then levels off at f1,
// rather than at t = 0.
// It needs Quadratic.
//
// Equivalent to:
//   - scipy.signal.chirp's vertex_zero=False
//   - MATLAB has none, as its other shape evaluates the sweep at fliplr(-t),
//     which depends on how t is laid out
func VertexAtT1() (option Option) {
	option = func(c *config) {
		c.vertexAtT1 = true
		c.given[optionVertexAtT1] = true
	}
	return
}

// Phase starts the waveform's cosine at phi degrees rather than 0, so 90 is a sine wave and -90 its negative.
//
// Units: phi is in degrees
//
// Equivalent to:
//   - scipy.signal.chirp's phi
//   - MATLAB chirp's phi
func Phase(phi float64) (option Option) {
	option = func(c *config) {
		c.phase = phi
		c.given[optionPhase] = true
	}
	return
}

// optionKind identifies an Option, for checking which apply to a waveform.
type optionKind int

const (
	optionLinear optionKind = iota
	optionQuadratic
	optionLogarithmic
	optionHyperbolic
	optionVertexAtT1
	optionPhase
	optionCount // the number of options, keep last
)

// optionNames names each option, for error messages.
var optionNames = [optionCount]string{
	optionLinear:      "Linear",
	optionQuadratic:   "Quadratic",
	optionLogarithmic: "Logarithmic",
	optionHyperbolic:  "Hyperbolic",
	optionVertexAtT1:  "VertexAtT1",
	optionPhase:       "Phase",
}

// optionSet is a set of options, where optionSet[o] is whether o is in it.
type optionSet [optionCount]bool

// without names the options in s that aren't in other.
func (s optionSet) without(other optionSet) (names []string) {
	for o := range optionCount {
		if s[o] && !other[o] {
			names = append(names, optionNames[o])
		}
	}
	return
}

// sweep is how Chirp's frequency changes from f0 to f1.
type sweep int

const (
	sweepLinear sweep = iota
	sweepQuadratic
	sweepLogarithmic
	sweepHyperbolic
)

// config is a waveform's settings, built from Options.
type config struct {
	sweep      sweep
	vertexAtT1 bool
	phase      float64   // in degrees
	given      optionSet // which options were passed, so a waveform can reject the rest
}

// newConfig builds a waveform's settings from opts,
// returning ErrInvalidOption for any option not in accepted.
func newConfig(waveform string, opts []Option, accepted optionSet) (cfg config, err error) {
	for _, opt := range opts {
		opt(&cfg)
	}
	if rejected := cfg.given.without(accepted); len(rejected) > 0 {
		err = fmt.Errorf("%w: %s doesn't apply to %s", ErrInvalidOption, strings.Join(rejected, ", "), waveform)
	}
	return
}
