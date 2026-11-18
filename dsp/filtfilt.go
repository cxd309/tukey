package dsp

import (
	"fmt"

	"github.com/cxd309/tukey/internal/floats"
)

// FiltFilt applies the filter b/a to x twice, forward then backward, for zero phase:
// features stay where they are in the signal,
// and the magnitude response is squared, so a -3 dB edge becomes -6 dB.
// x is extended past each end, by PaddingOdd and SciPy's padding length 3·max(len(a), len(b)) unless the options say otherwise,
// and each pass starts from the filter's steady state scaled to its first sample,
// so neither end shows a startup transient.
// The filter must be stable, which isn't checked, same as SciPy and MATLAB.
//
//	f, err := dsp.Butter(2, dsp.Lowpass(0.3))
//	b, a, err := f.BA()
//	y, err := dsp.FiltFilt(b, a, x)
//
// Units:
//   - b and a are dimensionless
//   - y is in x's unit
//
// Returns:
//   - ErrInvalidCoefficients for an empty b or a, an a[0] of 0, or a NaN or infinite coefficient.
//   - ErrInvalidParameter for a negative PadLen.
//   - ErrSignalTooShort for an x no longer than the padding length.
//   - ErrNoSteadyState for a filter with a pole at DC, where Σ a = 0.
//
// Equivalent to:
//   - scipy.signal.filtfilt(b, a, x)
//   - MATLAB filtfilt(b, a, x), with MATLABPadLen
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

// SOSFiltFilt applies the filter sos to x twice, forward then backward, for zero phase,
// staying accurate at high orders and narrow bands.
// The magnitude response is squared, so a -3 dB edge becomes -6 dB.
// x is extended past each end, by PaddingOdd and SciPy's padding length 3·ntaps unless the options say otherwise,
// where ntaps is 2·len(sos) + 1, less 1 for each first-order section, whose b2 and a2 are both 0.
// Each pass starts from the filter's steady state scaled to its first sample,
// so neither end shows a startup transient.
// By default the whole cascade is filtered forward then backward, same as SciPy,
// and with BySection each section is in turn, same as MATLAB:
//
//	y, err := dsp.SOSFiltFilt(sos, x, dsp.BySection(), dsp.MATLABPadLen()) // same as MATLAB filtfilt(sos, 1, x)
//
// Units:
//   - sos is dimensionless
//   - y is in x's unit
//
// Returns:
//   - ErrInvalidSOS for no sections, a section whose a0 isn't 1, or a NaN or infinite coefficient.
//   - ErrInvalidParameter for a negative PadLen.
//   - ErrSignalTooShort for an x no longer than the padding length.
//   - ErrNoSteadyState for a filter with a pole at DC, where Σ a = 0.
//
// Equivalent to:
//   - scipy.signal.sosfiltfilt(sos, x)
//   - MATLAB filtfilt(sos, g, x), with BySection and MATLABPadLen, and g folded into the first section
func SOSFiltFilt(sos SOS, x []float64, opts ...FiltOption) (y []float64, err error) {
	c, err := sos.toCascade()
	if err != nil {
		return nil, err
	}
	y, err = c.zeroPhase(x, sos.size(), opts)
	return
}

// FiltOption configures FiltFilt's and SOSFiltFilt's edge handling:
// PadLen, MATLABPadLen, BySection, and the paddings PaddingOdd, PaddingEven, PaddingConstant and PaddingNone.
type FiltOption func(*filtConfig)

// PadLen sets how many samples FiltFilt and SOSFiltFilt add beyond each end of x, which must be longer than n.
// The default is SciPy's 3·ntaps, and each function's doc says how it counts ntaps.
//
// Units: n is a number of samples
//
// Equivalent to:
//   - scipy.signal.filtfilt's and sosfiltfilt's padlen=n
//   - MATLAB has none, as its length is fixed
func PadLen(n int) (option FiltOption) {
	option = func(c *filtConfig) {
		c.padLen = func(filterSize) int { return n }
	}
	return
}

// MATLABPadLen uses MATLAB's padding length, 3 times the filter's order and at least 1,
// in place of SciPy's.
// The order counts b and a up to their last nonzero coefficients,
// and for SOSFiltFilt each section's, summed.
// So an odd-order Butterworth's first section, its real pole with two zeros, is order 2.
// Both forms are checked against MATLAB R2025a.
//
// Equivalent to:
//   - SciPy's padlen=3*order
//   - MATLAB filtfilt(b, a, x) with FiltFilt, and filtfilt(sos, g, x) with SOSFiltFilt and BySection
func MATLABPadLen() (option FiltOption) {
	option = func(c *filtConfig) {
		c.padLen = func(size filterSize) int { return max(1, 3*size.order) }
	}
	return
}

// BySection makes SOSFiltFilt filter each section forward and backward in turn,
// each padded by the whole filter's padding length,
// rather than the whole cascade forward then backward.
// The two differ by their edge transients, which for a slow filter reach well into the signal.
// A FiltFilt filter is one section, so BySection doesn't change it.
//
// Equivalent to:
//   - SciPy has none, as sosfiltfilt filters the whole cascade
//   - MATLAB filtfilt(sos, g, x), with MATLABPadLen
func BySection() (option FiltOption) {
	option = func(c *filtConfig) {
		c.bySection = true
	}
	return
}

// PaddingOdd extends x past each end by reflecting it through its end sample,
// so the slope continues across the end.
// It's the default.
//
// Equivalent to:
//   - scipy.signal.filtfilt's padtype='odd', its default
//   - MATLAB filtfilt, which has no other
func PaddingOdd() (option FiltOption) {
	option = func(c *filtConfig) { c.padding = paddingOdd }
	return
}

// PaddingEven extends x past each end by mirroring it about its end sample,
// so the slope flips across the end.
//
// Equivalent to:
//   - scipy.signal.filtfilt's padtype='even'
//   - MATLAB has none
func PaddingEven() (option FiltOption) {
	option = func(c *filtConfig) { c.padding = paddingEven }
	return
}

// PaddingConstant extends x past each end by repeating its end sample,
// so the slope drops to 0 across the end.
//
// Equivalent to:
//   - scipy.signal.filtfilt's padtype='constant'
//   - MATLAB has none
func PaddingConstant() (option FiltOption) {
	option = func(c *filtConfig) { c.padding = paddingConstant }
	return
}

// PaddingNone filters x as given, with no extension, so the startup transients land on its ends.
// The padding length is then ignored.
//
// Equivalent to:
//   - scipy.signal.filtfilt's padtype=None
//   - MATLAB has none
func PaddingNone() (option FiltOption) {
	option = func(c *filtConfig) { c.padding = paddingNone }
	return
}

// padding is how FiltFilt extends the signal past each end before filtering,
// so each pass's startup transient falls on the extension rather than the data.
type padding int

const (
	paddingOdd padding = iota
	paddingEven
	paddingConstant
	paddingNone
)

// String names the padding by its option, for error messages.
func (p padding) String() (name string) {
	switch p {
	case paddingOdd:
		name = "PaddingOdd"
	case paddingEven:
		name = "PaddingEven"
	case paddingConstant:
		name = "PaddingConstant"
	case paddingNone:
		name = "PaddingNone"
	default:
		name = fmt.Sprintf("padding(%d)", int(p))
	}
	return
}

// filtConfig is FiltFilt's edge handling, built from FiltOptions.
type filtConfig struct {
	padding padding
	// padLen gives the padding length for a filter of the given size,
	// and nil means SciPy's default, 3·ntaps
	padLen func(size filterSize) (n int)
	// bySection filters each section forward and backward in turn
	bySection bool
}

// filterSize is a filter's size, counted both ways a padding length is:
// ntaps, same as SciPy, for its default 3·ntaps,
// and order, same as MATLAB's filtfilt, for its 3·order.
type filterSize struct {
	ntaps int
	order int
}

// effectiveOrder is the order of the filter b/a as MATLAB's filtfilt counts it:
// the index of the later of b's and a's last nonzero coefficients, so trailing zeros don't count.
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

// padLength resolves the padding length for a filter of the given size, 3·ntaps by default, same as SciPy,
// and checks it against the signal's length.
func (cfg filtConfig) padLength(size filterSize, signalLen int) (n int, err error) {
	n = 3 * size.ntaps
	if cfg.padLen != nil {
		n = cfg.padLen(size)
	}
	if cfg.padding == paddingNone {
		n = 0 // no padding means padlen is ignored, same as SciPy
	}
	if n < 0 {
		return 0, fmt.Errorf("%w: padding length must be >= 0, got %d", ErrInvalidParameter, n)
	}
	if signalLen <= n {
		return 0, fmt.Errorf("%w: len(x) = %d must be greater than the padding length %d", ErrSignalTooShort, signalLen, n)
	}
	return
}

// zeroPhase filters x forward and backward through the cascade, padded as opts say,
// the pipeline shared by FiltFilt and SOSFiltFilt.
// size sets the padding length, which with BySection is the same for every section.
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
func (c cascade) forwardBackward(x []float64, n int, padding padding) (y []float64, err error) {
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

// extend returns a new slice holding x with n extension samples before and after it,
// same as scipy.signal's odd_ext, even_ext and const_ext.
// It panics unless 0 ≤ n < len(x), and n = 0 for PaddingNone.
// FiltFilt validates its input first, so a panic here means a bug inside dsp.
func (p padding) extend(x []float64, n int) (extended []float64) {
	if !(n >= 0 && n < len(x)) {
		panic(fmt.Sprintf("dsp: extend needs 0 <= n < len(x), got n=%d, len(x)=%d", n, len(x)))
	}
	if p == paddingNone && n != 0 {
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

// reflect is the extension's sample that mirrors the sample v about the end sample edge.
func (p padding) reflect(edge, v float64) (sample float64) {
	switch p {
	case paddingOdd:
		sample = 2*edge - v
	case paddingEven:
		sample = v
	case paddingConstant:
		sample = edge
	default:
		panic(fmt.Sprintf("dsp: unknown padding %v", p))
	}
	return
}
