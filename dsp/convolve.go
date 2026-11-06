package dsp

import (
	"fmt"
	"math"

	"github.com/cxd309/tukey/fft"
	"github.com/cxd309/tukey/internal/floats"
)

// ConvolveMethod is how Convolve and Correlate compute their sums.
type ConvolveMethod int

const (
	// ConvolveAuto chooses ConvolveDirect or ConvolveFFT by which is faster for the inputs' lengths.
	// The default option
	ConvolveAuto ConvolveMethod = iota
	// ConvolveDirect sums every product directly, each kept exactly and summed with compensation,
	// so each output is within about 1 ULP of its exact value, however small it is beside the others.
	// O(len(a)·len(v)).
	ConvolveDirect
	// ConvolveFFT multiplies the inputs' spectra, O(n log n).
	// Its error is relative to the largest output, so outputs much smaller than that lose relative accuracy:
	// 1e-3 noise beside a 1e8 spikes keeps only a few digits.
	ConvolveFFT
)

// ConvolveOption configures Convolve and Correlate
type ConvolveOption func(*convolveConfig)

// convolveOutput is which part of the full result is returned
type convolveOutput int

// bounds is where the selected outputs lie in the full result of na and nv samples
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

const (
	outputFull convolveOutput = iota
	outputSame
	outputValid
)

// convolveConfig is Convolve's and Correlate's settings, built from ConvolveOptions.
type convolveConfig struct {
	output convolveOutput
	method ConvolveMethod
}

// newConvolveConfig applies opts and checks the inputs, for the function name
func newConvolveConfig(name string, a, v []float64, opts []ConvolveOption) (cfg convolveConfig, err error) {
	for _, opt := range opts {
		opt(&cfg)
	}
	switch {
	case len(a) == 0 || len(v) == 0:
		err = fmt.Errorf("%w: %s of %d and %d samples, needs at least 1 each", ErrInvalidLength, name, len(a), len(v))
	case cfg.method < ConvolveAuto || cfg.method > ConvolveFFT:
		err = fmt.Errorf("%w: ConvolveMethod %d", ErrInvalidMode, cfg.method)
	}
	return
}

// ConvolveSame returns the central len(a) outputs of the full result,
// from (len(v)-1)/2 rounded down, same numpy and scipy "same" mode
func ConvolveSame() ConvolveOption {
	return func(c *convolveConfig) { c.output = outputSame }
}

// ConvolveValid returns only the outputs where the shorter input lies wholly inside the longer,
// |len(a) - len(v)| + 1 of them, same as scipy "valid" mode.
// Convolution is commutative, so it doesn't matter which input is longer.
func ConvolveValid() ConvolveOption {
	return func(c *convolveConfig) { c.output = outputValid }
}

// ConvolveBy sets how the sums are computed. ConvolveAuto by default.
func ConvolveBy(method ConvolveMethod) ConvolveOption {
	return func(c *convolveConfig) { c.method = method }
}

// Convolve returns the convolution of a and v, c[k] = Σ a[i]·v[k-i]:
// by default the full results, len(a) + len(v) - 1 outputs,
// or the segment selected by ConvolveSame or ConvolveValid.
// Filtering a signal by an FIR filter h is
//
//	Convolve(x, h)
//
// or to keep x's length
//
//	Convolve(x, h, ConvolveSame())
//
// Returns:
//   - ErrInvalidLength if a or v is empty
//   - ErrInvalidMode for an unknown ConvolveMethod
//
// Equivalent to:
//   - numpy.convolve(a, v, mode)
//   - scipy.signal.convolve(a, v, mode)
//   - matlab conv(a, v, shape), except that "same" starts one output later for an even len(v) and "valid" is empty when v is longer than a
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

// Correlation is a cross-correlation: Values[i] is at lag Lags[i]
type Correlation struct {
	Values []float64
	Lags   []int
}

// Correlate returns the cross-correlation of a and v, Σ a[n+k]·v[n] at each lag k:
// by default every lag at which they overlap, -(len(v)-1) to len(a)-1,
// or the part ConvolveSame or ConvolveValid selects.
// If a is v delayed by d samples, the largest value is at lag d.
//
// # Returns:
//   - ErrInvalidLength if a or v is empty
//   - ErrInvalidMode for an unknown ConvolveMethod
//
// # Equivalent to:
//   - scipy.signal.correlate(a, v, mode), with Lags from scipy.signal.correlation_lags
//
// except correlation_lags labels mode "same" one lag high for an odd len(a) and an even len(v)
//
// - MATLAB [r, lags] = xcorr(a, v) at the lags where a and v overlap
//
// xcorr pads the shorter with zeros and gives every lag up to the longer's length
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

// convolveFull is the full convolution of a and v, by method
func convolveFull(a, v []float64, method ConvolveMethod) (full []float64, err error) {
	if method == ConvolveAuto {
		method = fasterMethod(len(a), len(v))
	}
	if method == ConvolveFFT {
		return convolveFFT(a, v)
	}
	full = convolveDirect(a, v)
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

// convolveFFT multiplies a's and v's spectra,
// both zero-padded to at least the full result's length so the circular convolution doesn't wrap
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

// fftCostRatio is how many direct multiply-adds cost the same as one unit of an FFT's size·log2(size).
// Provisional, to be set by benchmark
const fftCostRatio = 20

// fasterMethod is ConvolveDirect or ConvolveFFT,
// whichever the cost model says is faster for inputs of na and nv samples.
func fasterMethod(na, nv int) (method ConvolveMethod) {
	size := float64(fft.FastLength(na + nv - 1))
	if float64(na)*float64(nv) <= fftCostRatio*size*math.Log2(size) {
		return ConvolveDirect
	}
	return ConvolveFFT
}
