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
//   - scipy.signal.chirp method='linear', its default
//   - MATLAB chirp "linear", its default
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
//   - scipy.signal.chirp method='quadratic'
//   - MATLAB chirp "quadratic", with its default shape, for the vertex at t = 0
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
//   - scipy.signal.chirp method='logarithmic'
//   - MATLAB chirp "logarithmic", except at f0 = f1, where MATLAB divides by 0
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
//   - scipy.signal.chirp method='hyperbolic'
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
//   - scipy.signal.chirp vertex_zero=False
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
//   - scipy.signal.chirp phi
//   - MATLAB chirp phi
func Phase(phi float64) (option Option) {
	option = func(c *config) {
		c.phase = phi
		c.given[optionPhase] = true
	}
	return
}

// Duty sets the fraction of each period a waveform spends in the first part of its cycle:
// Square's high part, 0.5 by default.
// It must be between 0 and 1, where 0 is never high and 1 always.
//
// Units: duty is a fraction of the period
//
// Equivalnet to:
//   - scipy.signal.square duty
//   - MATLAB square's duty, which is a percentage, so Duty(0.25) is MATLAB 25
func Duty(duty float64) (option Option) {
	option = func(c *config) {
		c.duty = duty
		c.given[optionDuty] = true
	}
	return
}

// Bandwidth sets GaussPulse's bandwidth as a fraction of its centre frequency fc:
// its spectrum falls to BandwidthLevel at fc·(1 ± bw/2).
// It's 0.5 by default, and must be above 0.
//
// Units: bw is a fraction of fc
//
// Equivalent to:
//   - scipy.signal.gausspulse bw
//   - MATLAB gauspuls bw
func Bandwidth(bw float64) (option Option) {
	option = func(c *config) {
		c.bandwidth = bw
		c.given[optionBandwidth] = true
	}
	return
}

// BandwidthLevel sets the level below its peak where GaussPulse's spectrum is Bandwidth wide.
// It's -6 dB by default, and must be below 0.
//
// Units: levelDB is in dB, negative
//
// Equivalent to:
//   - scipy.signal.gausspulse bwr
//   - MATLAB gauspuls bwr
func BandwidthLevel(levelDB float64) (option Option) {
	option = func(c *config) {
		c.bandwidthLevel = levelDB
		c.given[optionBandwidthLevel] = true
	}
	return
}

// CutoffLevel sets the level below its peak where GaussPulseCutoff finds the envelope.
// It's -60 dB by default, and must be below 0.
//
// Units: levelDB is in dB, negative
//
// Equivalent to:
//   - scipy.signal.gausspulse tpr
//   - MATLAB gauspuls tpr
func CutoffLevel(levelDB float64) (option Option) {
	option = func(c *config) {
		c.cutoffLevel = levelDB
		c.given[optionCutoffLevel] = true
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
	optionDuty
	optionBandwidth
	optionBandwidthLevel
	optionCutoffLevel
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
	optionDuty:        "Duty",
	optionBandwidth:   "BandwidthLevel",
	optionCutoffLevel: "CutoffLevel",
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
	sweep          sweep
	vertexAtT1     bool
	phase          float64   // in degrees
	duty           float64   // a fraction of the period, given only with duty
	bandwidth      float64   // a fraction of fc, given only with Bandwidth
	bandwidthLevel float64   // in dB, given only with BandwidthLevel
	cutoffLevel    float64   // in dB, given only with CutoffLevel
	given          optionSet // which options were passed, so a waveform can reject the rest
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
