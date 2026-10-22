package dsp

import (
	"fmt"
	"math"
	"math/cmplx"
	"slices"

	"github.com/cxd309/tukey/fft"
	"github.com/cxd309/tukey/internal/floats"
)

// SignalMatch is the window of data closest to the signal FindSignal looks for.
type SignalMatch struct {
	Start    int
	Stop     int // inclusive, same as MATLAB istop
	Distance float64
}

// SearchMethod is how FindSignal scores each window of data.
type SearchMethod int

const (
	// SearchHybrid scores every window by FFT correlation,
	// then rescores exactly each window that, within the correlation's rounding,
	// could be the closest.
	// Direct's answer at close to Fast's cost. The default search method.
	SearchHybrid SearchMethod = iota
	// SearchDirect scores every window exactly,
	// summing its squared differences with compensation. O(len(data) × len(signal)).
	// Simplest, but slowest at every size measured:
	// 77 s where Hybrid takes 29 ms, for a 30,000-sample signal in 400,000 samples.
	SearchDirect
	// SearchFast scores every window by FFT correlation alone,
	// which rounds the same way as MATLAB findsignal:
	// a perfect match may not score exactly 0,
	// and windows whose distances differ by less than that rounding may be confused.
	// 2 to 5% faster than Hybrid, and about 2x faster on near-periodic data,
	// where many windows are nearly as close as the best.
	SearchFast
)

// FindSignalOption configures FindSignal.
type FindSignalOption func(*findSignalConfig)

// findSignalConfig is FindSignal's settings, built from FindSignalOptions.
type findSignalConfig struct {
	method SearchMethod
}

// SignalSearch sets how FindSignal scores each window.
func SignalSearch(method SearchMethod) FindSignalOption {
	return func(c *findSignalConfig) { c.method = method }
}

// FindSignal returns the window of data closest to signal, by the sum of
// squared differences, with ties to the earliest window.
//
// Returns:
//   - ErrInvalidLength for an empty signal, or one longer than data,
//   - ErrInvalidMode for an unknown SearchMethod.
//
// Equivalent to MATLAB [istart, istop, dist] = findsignal(data, signal),
// with 0-based indices. No SciPy equivalent.
func FindSignal(data, signal []float64, opts ...FindSignalOption) (match SignalMatch, err error) {
	if len(signal) == 0 || len(signal) > len(data) {
		return match, fmt.Errorf("%w: signal of length %d in data of length %d", ErrInvalidLength, len(signal), len(data))
	}
	var cfg findSignalConfig
	for _, opt := range opts {
		opt(&cfg)
	}
	switch cfg.method {
	case SearchDirect:
		match = findSignalDirect(data, signal)
	case SearchFast:
		match, err = findSignalFast(data, signal)
	case SearchHybrid:
		match, err = findSignalHybrid(data, signal)
	default:
		err = fmt.Errorf("%w: SearchMethod %d", ErrInvalidMode, cfg.method)
	}
	return
}

// findSignalDirect scores every window of data exactly,
// keeping the first of any equally close.
func findSignalDirect(data, signal []float64) (match SignalMatch) {
	m := len(signal)
	match.Distance = windowDistance(data[:m], signal)
	for start := 1; start+m <= len(data); start++ {
		if d := windowDistance(data[start:start+m], signal); d < match.Distance {
			match.Start, match.Distance = start, d
		}
	}
	match.Stop = match.Start + m - 1
	return
}

// windowDistance is the sum of squared differences between window and
// signal, summed with compensation.
func windowDistance(window, signal []float64) (distance float64) {
	var sum floats.Sum
	for i, v := range window {
		diff := v - signal[i]
		sum.Add(diff * diff)
	}
	distance = sum.Total()
	return
}

// findSignalFast scores every window by distanceProfile,
// keeping the first of any equally close.
func findSignalFast(data, signal []float64) (match SignalMatch, err error) {
	profile, err := distanceProfile(data, signal)
	if err != nil {
		return match, err
	}
	for i, d := range profile {
		if d < profile[match.Start] {
			match.Start = i
		}
	}
	// a sum of squares can't be negative, though the expanded form can round below 0
	match.Distance = max(profile[match.Start], 0)
	match.Stop = match.Start + len(signal) - 1
	return
}

// distanceProfile is every window's distance from signal, by the expanded form
// sum d^2 - 2 sum dt + sum t^2, with the cross term from an FFT correlation.
// It rounds at the scale of the sums of squares rather than of each distance,
// so a perfect match can come out slightly above or below 0.
func distanceProfile(data, signal []float64) (profile []float64, err error) {
	m := len(signal)
	cross, err := crossCorrelation(data, signal)
	if err != nil {
		return nil, err
	}
	var signalSquares, windowSquares floats.Sum
	for _, v := range signal {
		signalSquares.Add(v * v)
	}
	for _, v := range data[:m] {
		windowSquares.Add(v * v)
	}
	profile = make([]float64, len(data)-m+1)
	for i := range profile {
		if i > 0 { // slide the window on one sample
			windowSquares.Add(-data[i-1] * data[i-1])
			windowSquares.Add(data[i+m-1] * data[i+m-1])
		}
		profile[i] = windowSquares.Total() - 2*cross[i] + signalSquares.Total()
	}
	return
}

// crossCorrelation is sum data[i+j]·signal[j] for every window start i, by FFT:
// data's spectrum times the conjugate of signal's. Padding both to at least
// len(data) keeps the circular correlation from wrapping into any window.
func crossCorrelation(data, signal []float64) (cross []float64, err error) {
	size := fft.FastLength(len(data))
	dataSpectrum, err := fft.ForwardReal(data, fft.N(size))
	if err != nil {
		return nil, err
	}
	signalSpectrum, err := fft.ForwardReal(signal, fft.N(size))
	if err != nil {
		return nil, err
	}
	for k := range dataSpectrum {
		dataSpectrum[k] *= cmplx.Conj(signalSpectrum[k])
	}
	full, err := fft.InverseReal(dataSpectrum, size)
	if err != nil {
		return nil, err
	}
	cross = full[:len(data)-len(signal)+1]
	return
}

// findSignalHybrid scores every window by distanceProfile, then rescores exactly
// each window that, within the profile's rounding, could be the closest. Its
// start and distance are then Direct's, at close to Fast's cost.
func findSignalHybrid(data, signal []float64) (match SignalMatch, err error) {
	profile, err := distanceProfile(data, signal)
	if err != nil {
		return match, err
	}
	cutoff := slices.Min(profile) + 2*profileErrorBound(data, signal)
	m := len(signal)
	match.Start = -1
	for i, p := range profile {
		if p > cutoff {
			continue
		}
		if d := windowDistance(data[i:i+m], signal); match.Start < 0 || closer(d, match.Distance) {
			match.Start, match.Distance = i, d
		}
	}
	match.Stop = match.Start + m - 1
	return
}

// profileErrorBound bounds distanceProfile's rounding error in any window.
// The FFT correlation's error scales with the whole of data and signal,
// and grows with the transform's depth.
// The final additions' error is within the same bound,
// as each window's sums of squares are within the whole's.
// The constant is generous: too loose only costs extra windows rescored,
// where too tight could miss the closest.
// Only the bound's size matters, so its sums are plain, not compensated.
// TestProfileErrorBound checks it holds.
func profileErrorBound(data, signal []float64) (bound float64) {
	var squares float64
	for _, v := range data {
		squares += v * v
	}
	for _, v := range signal {
		squares += v * v
	}
	depth := math.Log2(float64(fft.FastLength(len(data)))) + 1
	bound = 16 * floats.Epsilon * depth * squares
	return
}

// closer reports whether distance d beats best.
// A NaN distance never wins, and anything that isn't NaN beats a NaN best,
// so NaN in the data can't hide a later match.
func closer(d, best float64) (ok bool) {
	ok = d < best || (math.IsNaN(best) && !math.IsNaN(d))
	return
}
