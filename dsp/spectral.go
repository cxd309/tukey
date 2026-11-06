package dsp

import (
	"fmt"
	"math"
	"slices"
	"strings"

	"github.com/cxd309/tukey/fft"
	"github.com/cxd309/tukey/internal/floats"
)

// spectralOptionKind identifies a SpectralOption, to check which apply to a function
type spectralOptionKind int

const (
	optionOverlap spectralOptionKind = iota
	optionNFFT
	optionAtFrequencies
	optionTwoSided
	optionSegmentDetrend
	optionScaling
	optionAverage
	spectralOptionCount // the number of options, keep last
)

// SpectralScaling is the units power is given in
type SpectralScaling int

const (
	// ScalingDensity gives power spectral density, in units^2/Hz:
	// comparable across samples rates and windows, and the default
	ScalingDensity SpectralScaling = iota
	// ScalingSpectrum gives power spectrum, in units^2: a sine's peak reads its power
	ScalingSpectrum
)

// SpectralAverage is how Welch combines its segments' spectra
type SpectralAverage int

const (
	// AverageMean averages them, default option
	AverageMean SpectralAverage = iota
	// AverageMedian takes their median, corrected for its bias,
	// which is robust to a few segments with transients in them
	AverageMedian
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

// spectralOptionNames names each option, for error messages.
var spectralOptionNames = [spectralOptionCount]string{
	optionOverlap:        "Overlap",
	optionNFFT:           "NFFT",
	optionAtFrequencies:  "AtFrequencies",
	optionTwoSided:       "TwoSided",
	optionSegmentDetrend: "SegmentDetrend",
	optionScaling:        "Scaling",
	optionAverage:        "Average",
}

// spectralConfig is a spectral function's settings built from SpectralOptions
type spectralConfig struct {
	overlap     int
	nfft        int
	frequencies []float64
	twoSided    bool
	detrend     DetrendMode
	given       [spectralOptionCount]bool // which options were passed
	scaling     SpectralScaling
	average     SpectralAverage
}

// Spectrum is a power spectrum: one value per frequency
type Spectrum struct {
	Frequencies []float64 // Hz
	Power       []float64 // units^2/Hz, or units^2 with ScalingSpectrum
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
	tf.Frequencies = binFrequencies(cfg, fs)
	tf.Values, err = transformBins(windowed, cfg.nfft, cfg.twoSided)
	return

}

// Scaling sets the units power is given . ScalingDensity by default.
func Scaling(s SpectralScaling) SpectralOption {
	return func(c *spectralConfig) { c.scaling, c.given[optionScaling] = s, true }
}

// Average sets how Welch combines its segments' spectra. AverageMean by default
func Average(a SpectralAverage) SpectralOption {
	return func(c *spectralConfig) { c.average, c.given[optionAverage] = a, true }
}

// Spectrogram returns the power of x's STFT (see STFT): one power spectrum per segment
//
// Returns same as STFT. AtFrequencies doesn't apply
//
// Equivalent to scipy.signal.spectrogram(x, fs, window, detrend=False, mode='psd'), transposed.
// Equivalent to MATLAB [~, F, T, P] = spectrogram(x, window, nooverlap, nfft, fs), with P transposed and "power" for ScalingSpectrum
func Spectrogram(x []float64, fs float64, window []float64, opts ...SpectralOption) (tf TimeFrequency[float64], err error) {
	cfg, err := newSpectralConfig("Spectrogram", x, fs, window, opts, optionOverlap, optionNFFT, optionTwoSided, optionSegmentDetrend, optionScaling)
	if err != nil {
		return
	}
	power, times, err := segmentPowers(x, fs, window, cfg)
	if err != nil {
		return
	}
	tf = TimeFrequency[float64]{Frequencies: binFrequencies(cfg, fs), Times: times, Values: power}
	return
}

// Welch returns x's power spectrum by Welch's method:
// the average of its segments' power spectra, which is much less noisy than one periodogram,
// at the cost of the coarser frequency spacing of a segment's length.
//
// Returns as STFT does. AtFrequencies doesn't apply.
//
// Equivalent to scipy.signal.welch(x, fs, window, detrend=False).
// Equivalent with AverageMean to MATLAB [pxx, f] = pwelch(x, window, noverlap, nfft, fs), "power" for ScalingSpectrum
func Welch(x []float64, fs float64, window []float64, opts ...SpectralOption) (spectrum Spectrum, err error) {
	cfg, err := newSpectralConfig("Welch", x, fs, window, opts, optionOverlap, optionNFFT, optionTwoSided, optionSegmentDetrend, optionScaling, optionAverage)
	if err != nil {
		return
	}
	power, _, err := segmentPowers(x, fs, window, cfg)
	if err != nil {
		return
	}
	spectrum = Spectrum{Frequencies: binFrequencies(cfg, fs), Power: averaged(power, cfg.average)}
	return
}

// Periodogram returns x's power spectrum from one transform of the whole signal,
// multiplied by window, which must be as long as x.
// window.Boxcar(len(x)) gives the plain periodogram.
//
// Returns ErrInvalidLength for a window that isn't x's length, otherwise as STFT does.
// Overlap, AtFrequencies and Average don't apply.
//
// Equivalent to scipy.signal.periodogram(x, fs, window, detrend=False).
// Equivalent to MATLAB [pxx, f] = periodogram(x, window, nfft, fs), "power" for ScalingSpectrum.
func Periodogram(x []float64, fs float64, window []float64, opts ...SpectralOption) (spectrum Spectrum, err error) {
	if len(window) != len(x) {
		return spectrum, fmt.Errorf("%w: window of length %d for x of length %d, must be the same", ErrInvalidLength, len(window), len(x))
	}
	cfg, err := newSpectralConfig("Periodogram", x, fs, window, opts, optionNFFT, optionTwoSided, optionSegmentDetrend, optionScaling)
	if err != nil {
		return
	}
	power, _, err := segmentPowers(x, fs, window, cfg) // one segment: the whole signal
	if err != nil {
		return
	}
	spectrum = Spectrum{Frequencies: binFrequencies(cfg, fs), Power: power[0]}
	return
}

// binFrequencies is the frequency of each FFT bin, in Hz.
func binFrequencies(cfg spectralConfig, fs float64) (frequencies []float64) {
	if cfg.twoSided {
		return fft.Freq(cfg.nfft, 1/fs)
	}
	return fft.RealFreq(cfg.nfft, 1/fs)
}

// segmentPowers is each segment's power spectrum, scaled as cfg says.
// One-sided, each frequency is doubled to count the negative frequency it mirrors, except DC, and Nyquest when nfft is even, which have no mirror
func segmentPowers(x []float64, fs float64, window []float64, cfg spectralConfig) (power [][]float64, times []float64, err error) {
	windowed, times, err := segments(x, fs, window, cfg)
	if err != nil {
		return nil, nil, err
	}
	values, err := transformBins(windowed, cfg.nfft, cfg.twoSided)
	if err != nil {
		return nil, nil, err
	}
	var sum, squares floats.Sum
	for _, w := range window {
		sum.Add(w)
		squares.Add(w * w)
	}
	scale := 1 / (fs * squares.Total())
	if cfg.scaling == ScalingSpectrum {
		scale = 1 / (sum.Total() * sum.Total())
	}
	power = make([][]float64, len(values))
	for i, spectrum := range values {
		power[i] = make([]float64, len(spectrum))
		for k, v := range spectrum {
			p := (real(v)*real(v) + imag(v)*imag(v)) * scale
			if !cfg.twoSided && k > 0 && !(cfg.nfft%2 == 0 && k == cfg.nfft/2) {
				p *= 2
			}
			power[i][k] = p
		}
	}
	return
}

// averaged is the mean or median of each frequency's power across segments.
// The median is divided by its bias for exponentially distributed powers.
func averaged(power [][]float64, average SpectralAverage) (combined []float64) {
	combined = make([]float64, len(power[0]))
	column := make([]float64, len(power))
	for k := range combined {
		for i := range power {
			column[i] = power[i][k]
		}
		if average == AverageMedian {
			combined[k] = median(column) / medianBias(len(column))
			continue
		}
		var sum floats.Sum
		for _, p := range column {
			sum.Add(p)
		}
		combined[k] = sum.Total() / float64(len(column))
	}
	return
}

// median is the middle value of x, or the mean of the middle two. It sorts x
func median(x []float64) (m float64) {
	slices.Sort(x)
	n := len(x)
	if n%2 == 1 {
		return x[n/2]
	}
	return (x[n/2-1] + x[n/2]) / 2
}

// medianBias is how far below their mean the median of n exponentially distributed values falls:
// 1 - 1/2 + 1/3 - ... to the 1/n term for odd n, same as how scipy welch corrects
func medianBias(n int) (bias float64) {
	bias = 1
	for k := 1; k <= (n-1)/2; k++ {
		bias += 1/float64(2*k+1) - 1/float64(2*k)
	}
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
