package dsp

import (
	"fmt"
	"math"

	"github.com/cxd309/tukey/fft"
	"github.com/cxd309/tukey/internal/floats"
)

// Convolve returns the convolution of a and v, c[k] = Σ a[i]·v[k-i]:
// by default the full result, len(a) + len(v) - 1 outputs,
// or the part ConvolveSame or ConvolveValid selects.
// It sums directly or by FFT, whichever is faster, unless ConvolveDirect or ConvolveFFT says.
// Filtering a signal by an FIR filter h, keeping x's length, is:
//
//	y, err := dsp.Convolve(x, h, dsp.ConvolveSame())
//
// Units: c is in a's unit times v's
//
// Returns: ErrSignalTooShort for an empty a or v.
//
// Equivalent to:
//   - scipy.signal.convolve(a, v, mode), and numpy.convolve(a, v, mode)
//   - MATLAB conv(a, v, shape), except that its 'same' starts one output later for an even len(v),
//     and its 'valid' is empty when v is longer than a
func Convolve(a, v []float64, opts ...ConvolveOption) (c []float64, err error) {
	cfg, err := newConvolveConfig("Convolve", a, v, opts)
	if err != nil {
		return nil, err
	}
	full, err := convolveFull(a, v, cfg.method)
	if err != nil {
		return nil, err
	}
	start, stop := cfg.output.bounds(len(a), len(v))
	c = full[start:stop:stop]
	return
}

// Correlate returns the cross-correlation of a and v, Σ a[n+k]·v[n] at each lag k:
// by default every lag at which they overlap, -(len(v)-1) to len(a)-1,
// or the part ConvolveSame or ConvolveValid selects.
// If a is v delayed by d samples, the largest value is at lag d.
// It sums directly or by FFT, whichever is faster, unless ConvolveDirect or ConvolveFFT says.
//
//	corr, err := dsp.Correlate(received, sent)
//	delay := corr.Lags[slices.Index(corr.Values, slices.Max(corr.Values))]
//
// Units:
//   - Values is in a's unit times v's
//   - Lags are in samples
//
// Returns: ErrSignalTooShort for an empty a or v.
//
// Equivalent to:
//   - scipy.signal.correlate(a, v, mode), with Lags from scipy.signal.correlation_lags,
//     except that correlation_lags labels mode='same' one lag high for an odd len(a) and an even len(v)
//   - MATLAB [r, lags] = xcorr(a, v) at the lags where a and v overlap,
//     as xcorr pads the shorter with zeros and gives every lag up to the longer's length
func Correlate(a, v []float64, opts ...ConvolveOption) (c Correlation, err error) {
	cfg, err := newConvolveConfig("Correlate", a, v, opts)
	if err != nil {
		return c, err
	}
	full, err := convolveFull(a, floats.Reversed(v), cfg.method)
	if err != nil {
		return c, err
	}
	start, stop := cfg.output.bounds(len(a), len(v))
	c.Values = full[start:stop:stop]
	c.Lags = make([]int, stop-start)
	for i := range c.Lags {
		c.Lags[i] = start + i - (len(v) - 1) // full's first output is lag -(len(v)-1)
	}
	return
}

// ConvolveOption configures Convolve and Correlate:
// ConvolveSame or ConvolveValid for the output, and ConvolveAuto, ConvolveDirect or ConvolveFFT for the method.
type ConvolveOption func(*convolveConfig)

// ConvolveSame returns the central len(a) outputs of the full result, from (len(v)-1)/2 rounded down.
//
// Equivalent to:
//   - scipy.signal.convolve's mode='same', and numpy.convolve's
//   - MATLAB conv's 'same', except that it starts one output later for an even len(v)
func ConvolveSame() ConvolveOption {
	return func(c *convolveConfig) { c.output = outputSame }
}

// ConvolveValid returns only the outputs where the shorter input lies wholly inside the longer,
// |len(a) - len(v)| + 1 of them.
// Convolution is commutative, so it doesn't matter which input is longer.
//
// Equivalent to:
//   - scipy.signal.convolve's mode='valid'
//   - MATLAB conv's 'valid', except that it's empty when v is longer than a
func ConvolveValid() ConvolveOption {
	return func(c *convolveConfig) { c.output = outputValid }
}

// ConvolveAuto chooses ConvolveDirect or ConvolveFFT by which is faster for the inputs' lengths,
// which is ConvolveFFT for all but the shortest inputs.
// If either input has a NaN or infinite sample it uses ConvolveDirect,
// so that sample affects only the outputs it's summed into.
// It's the default.
// For inputs with a high dynamic range, use ConvolveDirect.
//
// Equivalent to:
//   - scipy.signal.convolve's method='auto', its default
//   - MATLAB has none, as conv always sums directly
func ConvolveAuto() ConvolveOption {
	return func(c *convolveConfig) { c.method = methodAuto }
}

// ConvolveDirect sums every product directly, each kept exactly and summed with compensation,
// so each output is within about 1 ULP of its exact value, however small it is beside the others.
// It takes O(len(a)·len(v)).
//
// Equivalent to:
//   - scipy.signal.convolve's method='direct', which sums plainly
//   - MATLAB conv, which sums plainly
func ConvolveDirect() ConvolveOption {
	return func(c *convolveConfig) { c.method = methodDirect }
}

// ConvolveFFT multiplies the inputs' spectra, in O(n log n).
// Its error is relative to the largest output, so outputs much smaller than that lose relative accuracy,
// e.g. 1e-3 noise beside 1e8 spikes keeps only a few digits.
// A NaN or infinite sample spreads to every output.
//
// Equivalent to:
//   - scipy.signal.convolve's method='fft', and scipy.signal.fftconvolve
//   - MATLAB has none
func ConvolveFFT() ConvolveOption {
	return func(c *convolveConfig) { c.method = methodFFT }
}

// Correlation is a cross-correlation, where Values[i] is at lag Lags[i].
type Correlation struct {
	Values []float64 // in a's unit times v's
	Lags   []int     // in samples, a positive lag meaning a is later than v
}

// convolveMethod is how Convolve and Correlate compute their sums.
type convolveMethod int

const (
	methodAuto convolveMethod = iota
	methodDirect
	methodFFT
)

// convolveConfig is Convolve's and Correlate's settings, built from ConvolveOptions.
type convolveConfig struct {
	output convolveOutput
	method convolveMethod
}

// convolveOutput is which part of the full result is returned.
type convolveOutput int

const (
	outputFull convolveOutput = iota
	outputSame
	outputValid
)

// bounds is where the selected outputs lie in the full result of na and nv samples.
func (o convolveOutput) bounds(na, nv int) (start, stop int) {
	switch o {
	case outputSame:
		start = (nv - 1) / 2
		stop = start + na
	case outputValid:
		start, stop = min(na, nv)-1, max(na, nv)
	default:
		start, stop = 0, na+nv-1
	}
	return
}

// newConvolveConfig applies opts and checks the inputs, returning ErrSignalTooShort for the function called name.
func newConvolveConfig(name string, a, v []float64, opts []ConvolveOption) (cfg convolveConfig, err error) {
	for _, opt := range opts {
		opt(&cfg)
	}
	switch {
	case len(a) == 0 || len(v) == 0:
		err = fmt.Errorf("%w: %s of %d and %d samples, needs at least 1 each", ErrSignalTooShort, name, len(a), len(v))
	}
	return
}

// convolveFull is the full convolution of a and v, by method.
func convolveFull(a, v []float64, method convolveMethod) (full []float64, err error) {
	if method == methodAuto {
		method = autoMethod(a, v)
	}
	if method == methodFFT {
		return convolveFFT(a, v)
	}
	full = convolveDirect(a, v)
	return
}

// autoMethod is the faster method for a and v,
// or ConvolveDirect if either has a NaN or infinite sample, which an FFT would spread to every output.
func autoMethod(a, v []float64) (method convolveMethod) {
	method = fasterMethod(len(a), len(v))
	if method == methodFFT && !(allFinite(a) && allFinite(v)) {
		method = methodDirect
	}
	return
}

// fasterMethod is ConvolveDirect or ConvolveFFT,
// whichever the cost model finds faster for inputs of na and nv samples.
func fasterMethod(na, nv int) (method convolveMethod) {
	size := float64(fft.FastLength(na + nv - 1))
	if float64(na)*float64(nv) <= fftCostOverhead+fftCostPerUnit*size*math.Log2(size) {
		return methodDirect
	}
	return methodFFT
}

// The FFT's measured cost, in units of one direct compensated product:
// a fixed overhead, then a cost per size·log2(size).
const (
	fftCostOverhead = 130
	fftCostPerUnit  = 0.4
)

// allFinite reports whether x has no NaN or infinite samples.
func allFinite(x []float64) bool {
	for _, s := range x {
		if math.IsNaN(s) || math.IsInf(s, 0) {
			return false
		}
	}
	return true
}

// convolveFFT multiplies a's and v's spectra,
// both zero-padded to at least the full result's length, so the circular convolution doesn't wrap.
func convolveFFT(a, v []float64) (full []float64, err error) {
	n := len(a) + len(v) - 1
	size := fft.FastLength(n)
	aSpectrum, err := fft.ForwardReal(a, fft.N(size))
	if err != nil {
		return nil, err
	}
	vSpectrum, err := fft.ForwardReal(v, fft.N(size))
	if err != nil {
		return nil, err
	}
	for k := range aSpectrum {
		aSpectrum[k] *= vSpectrum[k]
	}
	full, err = fft.InverseReal(aSpectrum, size)
	if err != nil {
		return nil, err
	}
	full = full[:n:n]
	return
}

// convolveDirect sums each output's products directly.
// Each product's rounding error is recovered exactly by FMA and summed with it,
// so each output is its exact sum of products, rounded close to once.
func convolveDirect(a, v []float64) (full []float64) {
	full = make([]float64, len(a)+len(v)-1)
	for k := range full {
		var sum floats.Sum
		for i := max(0, k-len(v)+1); i <= min(k, len(a)-1); i++ {
			p := float64(a[i] * v[k-i]) // the conversion rounds, so p can't be fused into the sum
			sum.Add(p)
			if !math.IsInf(p, 0) { // an infinite product has no finite rounding error
				sum.Add(math.FMA(a[i], v[k-i], -p))
			}
		}
		full[k] = sum.Total()
	}
	return
}
