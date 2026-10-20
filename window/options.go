package window

import (
	"fmt"
	"strings"
)

// optionKind identifies an Option, for checking which apply to a window
type optionKind int

const (
	optionPeriodic optionKind = iota
	optionCount               // the number of options, keep last
)

// optionNames names each option, for error messages
var optionNames = [optionCount]string{
	optionPeriodic: "Periodic",
}

// optionSet is a set of optionSet: optionSet[o] is whether o is in it
type optionSet [optionCount]bool

// Option configures a window (see Periodic)
// a window given an option that doesn't apply to it returns ErrInvalidParameter
type Option func(*config)

// config is a window's settings, built from Options
type config struct {
	periodic bool
	given    optionSet // which options were passed, so a window can reject the rest
}

// Periodic gives the DFT-even window, for spectral analysis:
// the symmetric window of length n+1 without its last point
// default is symmetric, for filter design
// equivalent to SciPy's sym=False and MATLAB's "periodic"
func Periodic() (option Option) {
	option = func(c *config) {
		c.periodic = true
		c.given[optionPeriodic] = true
	}
	return
}

// newConfig builds a window's condfig from opts,
// returning ErrInvalidParameter for any option not in accepted
func newConfig(window string, opts []Option, accepted optionSet) (cfg config, err error) {
	for _, opt := range opts {
		opt(&cfg)
	}
	if rejected := cfg.given.without(accepted); len(rejected) > 0 {
		err = fmt.Errorf("%w: %s doesn't apply to %s", ErrInvalidParameter, strings.Join(rejected, ", "), window)
	}
	return
}

// without names the options in s that aren't in other
func (s optionSet) without(other optionSet) (names []string) {
	for o := range optionCount {
		if s[0] && !other[o] {
			names = append(names, optionNames[o])
		}
	}
	return
}
