package dsp

import (
	"fmt"
)

// FiltOption configures FiltFilt's edge handling
// see PadLen, MATLABPadLen and PadType
type FiltOption func(*filtConfig)

// filtConfig is FiltFilt's edge handling, built up from FiltOptions
type filtConfig struct {
	padding Padding
	// padLen gives the padding length for a filter with ntaps coefficients
	// nil means SciPy's default, 3 * ntaps
	padLen func(ntaps int) (n int)
}

// PadLen sets how many samples FiltFilt adds beyond each end of x
// the default is SciPy's 3 * max(lan(a), len(b)), and x must be longer than n
func PadLen(n int) (option FiltOption) {
	option = func(c *filtConfig) {
		c.padLen = func(int) int { return n }
	}
	return
}

// MATLABPadLen uses the padding length MATLAB's filtfilt uses
// 3 * (max(len(a), len(b)) - 1), in place of SciPy's default
// only samples near the ends of the output are affected
func MATLABPadLen() (option FiltOption) {
	option = func(c *filtConfig) {
		c.padLen = func(ntaps int) int { return 3 * (ntaps - 1) }
	}
	return
}

// PadType sets how FiltFilt extends x past each end
// default is PaddingOdd
func PadType(p Padding) (option FiltOption) {
	option = func(c *filtConfig) {
		c.padding = p
	}
	return
}

// FiltFilt applies the IIR filter b, a to x twice, forward then backwards
// giving zeros phase distortion: features stay where they are in the signal
// and the magnitude response is squared (so the -3dB cutoff becomes -6dB)
//
// x is extended past each end (see PadType, PadLen) and each pass starts from
// the filter's steady state scaled to the edge sample, so neither end shows a
// startup transient
// x must be longer than the padding length
//
// a must describe a stable filter (all poles strictly inside the unit circle):
// this isn't checked, matching SciPy and MATLAB
//
//	b, a, err := dsp.Butter(2, dsp.Lowpass(0.3))
//	y, err := dsp.FiltFilt(b, a, x)
//
// returns ErrInvalidCoefficients if b or a is empty or a[0] is zero,
// ErrInvalidPadding for an unknown Padding or a negative padding length,
// ErrSignalTooShort if len(x) isn't greater than the padding length,
// and ErrNoSteadyState if the filter has a pole at DC (sum(a) == 0)
//
// equivalent to scipy.signal.filtfilt(b, a, x)
// or to MATLAB filtfilt(b, a, x) with MATLABPadLen
func FiltFilt(b, a, x []float64, opts ...FiltOption) (y []float64, err error) {
	f, err := newDigitalFilter(b, a)
	if err != nil {
		return nil, err
	}
	var cfg filtConfig
	for _, opt := range opts {
		opt(&cfg)
	}
	if !(cfg.padding >= PaddingOdd && cfg.padding <= PaddingNone) {
		return nil, fmt.Errorf("%w: unknown padding %v", ErrInvalidPadding, cfg.padding)
	}

	n := 3 * len(f.a) // SciPy's default, len(f.a) is max(len(a), len(b)) after newDigitalFilter
	if cfg.padLen != nil {
		n = cfg.padLen(len(f.a))
	}
	if cfg.padding == PaddingNone {
		n = 0 // as in SciPy, no padding means padlen is ignored
	}
	if n < 0 {
		return nil, fmt.Errorf("%w: padding length must be >= 0, got %d", ErrInvalidPadding, n)
	}
	if len(x) <= n {
		return nil, fmt.Errorf("%w: len(x) = %d must be greater than the padding length %d", ErrSignalTooShort, len(x), n)
	}

	zi, err := f.steadyState()
	if err != nil {
		return nil, err
	}

	extended := cfg.padding.extend(x, n)
	forward := f.apply(extended, scaled(zi, extended[0]))
	backward := f.apply(reversed(forward), scaled(zi, forward[len(forward)-1]))
	y = reversed(backward)[n : len(backward)-n]
	return
}
