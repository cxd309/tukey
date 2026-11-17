package window

import (
	"fmt"
	"strings"
)

// Option configures a window: Periodic, Center or Unnormalised.
// A window given an option that doesn't apply to it returns ErrInvalidOption.
type Option func(*config)

// Periodic gives the DFT-even window, for spectral analysis:
// the symmetric window of length n+1 without its last point.
// Without it a window is symmetric, for filter design.
// Every window but KaiserBesselDerived takes it.
//
// Equivalent to:
//   - SciPy's sym=False
//   - MATLAB's 'periodic', for the windows that have it
func Periodic() (option Option) {
	option = func(c *config) {
		c.periodic = true
		c.given[optionPeriodic] = true
	}
	return
}

// Center moves Exponential's peak to sample c, which may be fractional.
// It needs Periodic, as a symmetric window is centred.
//
// Units: c is a sample index
//
// Equivalent to:
//   - scipy.signal.windows.exponential's center, with sym=False
//   - MATLAB has none
func Center(c float64) (option Option) {
	option = func(cfg *config) {
		cfg.center = c
		cfg.given[optionCenter] = true
	}
	return
}

// Unnormalised leaves Taylor's peak where its formula puts it,
// rather than scaling it to 1 at the centre.
//
// Equivalent to:
//   - scipy.signal.windows.taylor's norm=False
//   - MATLAB taylorwin, which doesn't normalise
func Unnormalised() (option Option) {
	option = func(cfg *config) {
		cfg.unnormalised = true
		cfg.given[optionUnnormalised] = true
	}
	return
}

// optionKind identifies an Option, for checking which apply to a window.
type optionKind int

const (
	optionPeriodic optionKind = iota
	optionCenter
	optionUnnormalised
	optionCount // the number of options, keep last
)

// optionNames names each option, for error messages.
var optionNames = [optionCount]string{
	optionPeriodic:     "Periodic",
	optionCenter:       "Center",
	optionUnnormalised: "Unnormalised",
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

// config is a window's settings, built from Options.
type config struct {
	periodic     bool
	center       float64
	unnormalised bool
	given        optionSet // which options were passed, so a window can reject the rest
}

// newConfig builds a window's settings from opts,
// returning ErrInvalidOption for any option not in accepted.
func newConfig(window string, opts []Option, accepted optionSet) (cfg config, err error) {
	for _, opt := range opts {
		opt(&cfg)
	}
	if rejected := cfg.given.without(accepted); len(rejected) > 0 {
		err = fmt.Errorf("%w: %s doesn't apply to %s", ErrInvalidOption, strings.Join(rejected, ", "), window)
	}
	return
}
