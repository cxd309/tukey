package dsp

import (
	"fmt"
	"math"
	"strings"

	"github.com/cxd309/tukey/fft"
	"github.com/cxd309/tukey/internal/floats"
)

// TimeFrequency is a signal's spectrum over time: one spectrum per segment
type TimeFrequency[T complex128 | float64] struct {
	Frequencies []float64 // Hz, one for each value in a segment's spectrum
	Times       []float64 // s, the centre of each segment
	Values      [][]T     // Values[segment][frequency]
}

// SpectralOption configures STFT, Spectrogram, Welch and Periodogram.
// A function given an option that doesn't apply to it returns ErrInvalidParameter
type SpectralOption func(*spectralConfig)

// spectralOptionKing identifies a SpectralOption, to check which apply to a function
type spectralOptionKind int

const (
	optionOverlap spectralOptionKind = iota
	optionNFFT
	optionAtFrequencies
	optionTwoSided
	optionSegmentDetrend
	spectralOptionCount // the number of options, keep last
)

// spectralOptionNames names each option, for error messages.
var spectralOptionNames = [spectralOptionCount]string{
	optionOverlap:        "Overlap",
	optionNFFT:           "NFFT",
	optionAtFrequencies:  "AtFrequencies",
	optionTwoSided:       "TwoSided",
	optionSegmentDetrend: "SegmentDetrend",
}

// spectralConfig is a spectral function's settings built from SpectralOptions
type spectralConfig struct {
	overlap     int
	nfft        int
	frequencies []float64
	twoSided    bool
	detrend     DetrendMode
	given       [spectralOptionCount]bool // which options were passed
}

// Overlap sets how many samples neighbouring segments share.
// Half the window by default
func Overlap(n int) SpectralOption {
	return func(c *spectralConfig) { c.overlap, c.given[optionOverlap] = n, true }
}

// NFFT zero-pads each segment to n samples before transforming it,
// for more finely spaced frequencies. The window's length by default.
func NFFT(n int) SpectralOption {
	return func(c *spectralConfig) { c.nfft, c.given[optionNFFT] = n, true }
}

// AtFrequencies evaluates each segment's transform at these frequencies, in Hz,
// rather than at the FFT's bins: to look closely at a band without a huge NFFT.
// It can't be combined with NFFT or TwoSided.
func AtFrequencies(f []float64) SpectralOption {
	f = append([]float64(nil), f...)
	return func(c *spectralConfig) { c.frequencies, c.given[optionAtFrequencies] = f, true }
}

// TwoSided keeps the negative frequencies too, ordered as fft.Freq orders them.
// One-sided by default, as a real signal's negative frequencies mirror its positive ones.
func TwoSided() SpectralOption {
	return func(c *spectralConfig) { c.twoSided, c.given[optionTwoSided] = true, true }
}

// SegmentDetrend detrends each segment, before it's windowed, same as Detrend.
// No detrending by default: to remove the whole signal's mean or trend, call Detrend first.
func SegmentDetrend(mode DetrendMode) SpectralOption {
	return func(c *spectralConfig) { c.detrend, c.given[optionSegmentDetrend] = mode, true }
}

// STFT returns the short-time Fourier transofrm of x, sampled at fs Hz:
// x cut into segments of len(window) samples,
// each multiplied by window and tranformed.
// Each segment's tranform is plain, sum w[n]·x[n]·e^(-j2πfn/fs),
// with its phase measured from the segment's start.
// Divide by the window's sum for amplitudes.
// A prtial last segment is dropped.
//
// Returns:
//   - ErrInvalidLength for an empty window, one longer than x,
//     an NFFT shorter than the window, or AtFrequencies with no frequencies
//   - ErrInvalidParameter for an fs that isn't positive, an Overlap that's negative
//     or at least the window's length, or an option that doesn't apply
//
// Equivalent to MATLAB [S, F, T] = spectrogram(x, window, nooverlap, nfft, fs), with S transposed,
// Equivalent to scipy.signal.stft(x, fs, window, boundary=None, padded=False, detrend=False)
// multiplied by the window's sum.
func STFT(x []float64, fs float64, window []float64, opts ...SpectralOption) (tf TimeFrequency[complex128], err error) {
	cfg, err := newSpectralConfig("STFT", x, fs, window, opts, optionOverlap, optionNFFT, optionAtFrequencies, optionTwoSided, optionSegmentDetrend)
	if err != nil {
		return
	}
	windowed, times, err := segments(x, fs, window, cfg)
	if err != nil {
		return
	}
	tf.Times = times
	if cfg.given[optionAtFrequencies] {
		tf.Frequencies = append([]float64(nil), cfg.frequencies...)
		tf.Values = transformAt(windowed, cfg.frequencies, fs)
		return
	}
	if cfg.twoSided {
		tf.Frequencies = fft.Freq(cfg.nfft, 1/fs)
	} else {
		tf.Frequencies = fft.RealFreq(cfg.nfft, 1/fs)
	}
	tf.Values, err = transformBins(windowed, cfg.nfft, cfg.twoSided)
	return

}

// newSpectralConfig builds a spectral function's settings from opts,
// with the defaults filled in, and checks them and its inputs.
// accepted are the options that apply to funciton.
func newSpectralConfig(function string, x []float64, fs float64, window []float64, opts []SpectralOption, accepted ...spectralOptionKind) (cfg spectralConfig, err error) {
	cfg = spectralConfig{overlap: len(window) / 2, nfft: len(window)}
	for _, opt := range opts {
		opt(&cfg)
	}
	var accepts [spectralOptionCount]bool
	for _, o := range accepted {
		accepts[o] = true
	}
	var rejected []string
	for o := range spectralOptionCount {
		if cfg.given[o] && !accepts[o] {
			rejected = append(rejected, spectralOptionNames[o])
		}
	}
	switch {
	case len(rejected) > 0:
		err = fmt.Errorf("%w: %s doesn't apply to %s", ErrInvalidParameter, strings.Join(rejected, ", "), function)
	case !(fs > 0):
		err = fmt.Errorf("%w: fs %g, must be positive", ErrInvalidParameter, fs)
	case len(window) == 0 || len(window) > len(x):
		err = fmt.Errorf("%w: window of length %d for x of length %d", ErrInvalidLength, len(window), len(x))
	case cfg.overlap < 0 || cfg.overlap >= len(window):
		err = fmt.Errorf("%w: Overlap %d, must be from 0 to %d", ErrInvalidParameter, cfg.overlap, len(window)-1)
	case cfg.nfft < len(window):
		err = fmt.Errorf("%w: NFFT %d, shorter than the window's %d", ErrInvalidLength, cfg.nfft, len(window))
	case cfg.given[optionAtFrequencies] && (cfg.given[optionNFFT] || cfg.given[optionTwoSided]):
		err = fmt.Errorf("%w: AtFrequencies can't be combined with NFFT or TwoSided", ErrInvalidParameter)
	case cfg.given[optionAtFrequencies] && len(cfg.frequencies) == 0:
		err = fmt.Errorf("%w: AtFrequencies with no frequencies", ErrInvalidLength)
	}
	return
}

// segments cuts x into every whole segment of len(window) samples,
// cfg.overlap apart, detrends each if cfg says to, and applies the window.
// A partial last segment is dropped.
// times is each segment's centre, in seconds.
func segments(x []float64, fs float64, window []float64, cfg spectralConfig) (windowed [][]float64, times []float64, err error) {
	length := len(window)
	for start := 0; start+length <= len(x); start += length - cfg.overlap {
		segment := x[start : start+length]
		if cfg.given[optionSegmentDetrend] {
			if segment, err = Detrend(segment, cfg.detrend); err != nil {
				return nil, nil, err
			}
		}
		w := make([]float64, length)
		for i, v := range segment {
			w[i] = v * window[i]
		}
		windowed = append(windowed, w)
		times = append(times, (float64(start)+float64(length)/2)/fs)
	}
	return
}

// tranformBins is each windowed segment's FFT, zero-padded to nfft:
// the non-negative frequencies, or every frequency if twoSided
func transformBins(windowed [][]float64, nfft int, twoSided bool) (values [][]complex128, err error) {
	values = make([][]complex128, len(windowed))
	for i, segment := range windowed {
		if twoSided {
			values[i], err = fft.Forward(complexified(segment), fft.N(nfft))
		} else {
			values[i], err = fft.ForwardReal(segment, fft.N(nfft))
		}
		if err != nil {
			return nil, err
		}
	}
	return
}

// complexified is x as complex values
func complexified(x []float64) (z []complex128) {
	z = make([]complex128, len(x))
	for i, v := range x {
		z[i] = complex(v, 0)
	}
	return
}

// transformAt is each windowed segment's tranform at frequencies, in Hz:
// sum segment[n]·e^(-j2pifn/fs), each part summed with compensation.
func transformAt(windowed [][]float64, frequencies []float64, fs float64) (values [][]complex128) {
	length := len(windowed[0])
	cos := make([][]float64, len(frequencies))
	sin := make([][]float64, len(frequencies))
	for k, f := range frequencies { // the same for every segment
		cos[k], sin[k] = twiddles(f, fs, length)
	}
	values = make([][]complex128, len(windowed))
	for i, segment := range windowed {
		values[i] = make([]complex128, len(frequencies))
		for k := range frequencies {
			var re, im floats.Sum
			for n, v := range segment {
				re.Add(v * cos[k][n])
				im.Add(-v * sin[k][n])
			}
			values[i][k] = complex(re.Total(), im.Total())
		}
	}
	return
}

// twiddles is the cosine and sine of 2 pi fn / fs for n from 0 to length -1.
// Each angle is reduced to a fraction of a turn before it's rounded:
// f n is split into its rounded product and that product's rounding error,
// found exactly with FMA, and the product is reduced modulo fs,
// which only subtracts multiples of fs and so rounds nothing.
// Rounding 2 pi fn / fs directly,
// an angle of up to hundreds of radians, would cost about 1e-13
func twiddles(f, fs float64, length int) (cos, sin []float64) {
	cos = make([]float64, length)
	sin = make([]float64, length)
	for n := range length {
		product := f * float64(n)
		lost := math.FMA(f, float64(n), -product) // exactly fn - product
		turns := (math.Mod(product, fs) + lost) / fs
		sin[n], cos[n] = math.Sincos(2 * math.Pi * turns)
	}
	return
}
