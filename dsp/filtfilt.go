package dsp

import (
	"fmt"

	"github.com/cxd309/tukey/internal/floats"
)

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
// ErrInvalidOption for an unknown Padding, ErrInvalidParameter for a negative padding length,
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
	// len(f.a) is max(len(a), len(b)) after newDigitalFilter
	size := filterSize{ntaps: len(f.a), order: effectiveOrder(f.b, f.a)}
	y, err = cascade{f}.zeroPhase(x, size, opts)
	return
}

// SOSFiltFilt is FiltFilt for a filter given as second-order sections:
// it filters x forward then backward for zero phase, with the same edge handling options.
// The default padding length is scipy 3* ntaps, where ntaps is 2* len(sos)+ 1,
// minus one for each first-order section (b2 and a2 both zero).
//
// By default the whole cascade is filtered forward then backward, as scipy sosfiltfilt.
// With BySection each section is, in turn, same as matlab filtfilt(sos, g, x),
// which with MATLABPadLen it matches:
//
//	y, err := dsp.SOSFiltFilt(sos, x, dsp.BySection(), dsp.MATLABPadLen())
//
// matlab g is in the first section already.
//
// Returns:
//   - ErrInvalidSOS same as SOSFilter (if sos has no sections or any section's a0 isn't 1)
//   - ErrInvalidOption, ErrInvalidParameter, ErrSignalTooShort and ErrNoSteadyState same as FiltFilt
//
// Equivalent to:
//   - scipy.signal.sosfiltfilt(sos, x)
//   - matlab filtfilt(sos, 1, x), with BySection and MATLABPadLen
func SOSFiltFilt(sos SOS, x []float64, opts ...FiltOption) (y []float64, err error) {
	c, err := sos.toCascade()
	if err != nil {
		return nil, err
	}
	y, err = c.zeroPhase(x, sos.size(), opts)
	return
}

// FiltOption configures FiltFilt and SOSFiltFilt edge handling
// see PadLen, MATLABPadLen, PadType and BySection
type FiltOption func(*filtConfig)

// PadLen sets how many samples FiltFilt and SOSFiltFilt add beyond each end of x
// x must be longer than n
// the default is SciPy's 3*ntaps: see each function for how it counts ntaps
func PadLen(n int) (option FiltOption) {
	option = func(c *filtConfig) {
		c.padLen = func(filterSize) int { return n }
	}
	return
}

// MATLABPadLen uses the padding length matlab filtfilt uses,
// 3* the filter's order, and at least 1, in place of scipy default.
// The order counts b and a up to their last nonzero coefficients,
// and for SOSFiltFilt each section, summed:
// an odd-order Butterworth's first section, its real pole with two zeros, is order 2
//
// With FiltFilt this matches matlab filtfilt(b, a, x),
// and with SOSFiltFilt and BySection its filtfilt(sos, g, x),
// both verified against matlab R2025a
func MATLABPadLen() (option FiltOption) {
	option = func(c *filtConfig) {
		c.padLen = func(size filterSize) int { return max(1, 3*size.order) }
	}
	return
}

// BySection makes SOSFiltFilt each section forward and backward in turn,
// each padded by the whole filter's padding length, same as matlab filtfilt(sos, g, x),
// rather than the whole cascade forward then backward, same as scipy sosfiltfilt.
// The two differ by their edge transients,
// which for a slow filter reach well into the signal.
// With MATLABPadLen it matches matlab
//
// A FiltFilt filter is one section, so BySection doesn't change it.
func BySection() (option FiltOption) {
	option = func(c *filtConfig) {
		c.bySection = true
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

// Padding selects how FiltFilt extends the signal past each end before filtering
// so that each pass's startup transient falls on the extension rather than the data
type Padding int

const (
	// PaddingOdd reflect the signal through its end sample
	// so both value and slope continue smoothly
	// default option like SciPy and MATLAB
	PaddingOdd Padding = iota
	// PaddingEven mirror the signal about its end sample
	PaddingEven
	// PaddingConstant repeats the end sample
	PaddingConstant
	// PaddingNone filters the signal as given
	PaddingNone
)

func (p Padding) String() (name string) {
	switch p {
	case PaddingOdd:
		name = "PaddingOdd"
	case PaddingEven:
		name = "PaddingEven"
	case PaddingConstant:
		name = "PaddingConstant"
	case PaddingNone:
		name = "PaddingNone"
	default:
		name = fmt.Sprintf("Padding(%d)", int(p))
	}
	return
}

// filtConfig is FiltFilt's edge handling, built up from FiltOptions
type filtConfig struct {
	padding Padding
	// padLen gives the padding length for a filter with ntaps coefficients
	// nil means SciPy's default, 3 * ntaps
	padLen func(ntaps filterSize) (n int)
	// bySection filters each section forward and backward in turn
	bySection bool
}

// filterSize is a filter's size, counted in both the ways padding length is:
// ntaps same as scipy, for its default 3* ntaps,
// and order same as matlab filtfilt, for its 3* order
type filterSize struct {
	ntaps int
	order int
}

// effectiveOrder is the order of the filter b/a as matlab filtfilt counts it:
// the index of the later of b and a last nonzero coefficients,
// so trailing zeros don't count
func effectiveOrder(b, a []float64) (order int) {
	for _, coeffs := range [][]float64{b, a} {
		for i := len(coeffs) - 1; i > order; i-- {
			if coeffs[i] != 0 {
				order = i
				break
			}
		}
	}
	return
}

// padLength resolves the padding length for a filter with ntaps coefficients
// (3*ntaps by default, as SciPy) and checks it against the signal's length
func (cfg filtConfig) padLength(size filterSize, signalLen int) (n int, err error) {
	if !(cfg.padding >= PaddingOdd && cfg.padding <= PaddingNone) {
		return 0, fmt.Errorf("%w: unknown padding %v", ErrInvalidOption, cfg.padding)
	}
	n = 3 * size.ntaps
	if cfg.padLen != nil {
		n = cfg.padLen(size)
	}
	if cfg.padding == PaddingNone {
		n = 0 // as in SciPy, no padding means padlen is ignored
	}
	if n < 0 {
		return 0, fmt.Errorf("%w: padding length must be >= 0, got %d", ErrInvalidParameter, n)
	}
	if signalLen <= n {
		return 0, fmt.Errorf("%w: len(x) = %d must be greater than the padding length %d", ErrSignalTooShort, signalLen, n)
	}
	return
}

// zeroPhase filters x forward and backward through the cascade, padded per opts,
// the pipeline shared by FiltFilt and SOSFiltFilt.
// size sets the padding length, which with BySection is the same for every section
func (c cascade) zeroPhase(x []float64, size filterSize, opts []FiltOption) (y []float64, err error) {
	var cfg filtConfig
	for _, opt := range opts {
		opt(&cfg)
	}
	n, err := cfg.padLength(size, len(x))
	if err != nil {
		return nil, err
	}
	if !cfg.bySection {
		return c.forwardBackward(x, n, cfg.padding)
	}
	y = x
	for _, section := range c {
		if y, err = (cascade{section}).forwardBackward(y, n, cfg.padding); err != nil {
			return nil, err
		}
	}
	return
}

// forwardBackward runs the cascade forward over x padded by n samples, then backward over the result,
// each pass starting from the steady state scaled to its first sample, and trims the padding.
func (c cascade) forwardBackward(x []float64, n int, padding Padding) (y []float64, err error) {
	zi, err := c.steadyState()
	if err != nil {
		return nil, err
	}

	extended := padding.extend(x, n)
	forward, _ := c.apply(extended, scaledState(zi, extended[0]))
	backward, _ := c.apply(floats.Reversed(forward), scaledState(zi, forward[len(forward)-1]))
	y = floats.Reversed(backward)[n : len(backward)-n]
	return
}

// extend returns a new slice holding x with n extension samples before and after it
// it panics unless 0 <= n < len(x), and n == 0 for PaddingNone
// FiltFilt validates caller input first, so a panic here means a bug inside dsp
// equvalent to scipy.signal's odd_ext, even_ext and const_ext
func (p Padding) extend(x []float64, n int) (extended []float64) {
	if !(n >= 0 && n < len(x)) {
		panic(fmt.Sprintf("dsp: extend needs 0 <= n < len(x), got n=%d, len(x)=%d", n, len(x)))
	}
	if p == PaddingNone && n != 0 {
		panic(fmt.Sprintf("dsp: extend with PaddingNone needs n = 0, got n=%d", n))
	}
	first, last := x[0], x[len(x)-1]
	extended = make([]float64, 0, len(x)+2*n)
	for i := n; i >= 1; i-- {
		extended = append(extended, p.reflect(first, x[i]))
	}
	extended = append(extended, x...)
	for i := 1; i <= n; i++ {
		extended = append(extended, p.reflect(last, x[len(x)-1-i]))
	}
	return
}

// reflect is the extensions sample that mirrors sample v about the edge sample edge
func (p Padding) reflect(edge, v float64) (sample float64) {
	switch p {
	case PaddingOdd:
		sample = 2*edge - v
	case PaddingEven:
		sample = v
	case PaddingConstant:
		sample = edge
	default:
		panic(fmt.Sprintf("dsp: unknown padding %v", p))
	}
	return
}
