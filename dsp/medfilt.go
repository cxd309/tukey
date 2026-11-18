package dsp

import (
	"fmt"
	"math"
	"slices"

	"github.com/cxd309/tukey/internal/floats"
)

// MedFilt replaces each sample of x with the median of the n samples centred on it,
// which removes spikes without blurring steps.
// n must be odd.
// A window holding a NaN gives NaN, and an infinite sample is outvoted, as any outlier is.
// A window with an even count, which only MedFiltTruncate gives near the ends,
// takes the midpoint of its two middle values, rounded once, so y is exact.
// It keeps the window sorted as it slides, O(n) per sample.
//
//	y, err := dsp.MedFilt(x, 5)
//
// Units:
//   - n is a number of samples
//   - y is in x's unit
//
// Returns: ErrInvalidLength for an n that isn't odd and positive.
//
// Equivalent to:
//   - scipy.signal.medfilt(x, n), with MedFiltZeroPad, except at a NaN, where its result depends on how NaN sorts
//   - MATLAB medfilt1(x, n, 'truncate'), or medfilt1(x, n) with MedFiltZeroPad
func MedFilt(x []float64, n int, opts ...MedFiltOption) (y []float64, err error) {
	var cfg medFiltConfig
	for _, opt := range opts {
		opt(&cfg)
	}
	switch {
	case n < 1 || n%2 == 0:
		return nil, fmt.Errorf("%w: MedFilt needs an odd, positive window length, got %d", ErrInvalidLength, n)
	}

	half := n / 2
	window := sortedWindow{sorted: make([]float64, 0, min(n, len(x)))}
	for _, v := range x[:min(half+1, len(x))] {
		window.add(v)
	}
	y = make([]float64, len(x))
	for i := range x {
		if i > half {
			window.remove(x[i-half-1])
		}
		if next := i + half; i > 0 && next < len(x) {
			window.add(x[next])
		}
		zeros := 0
		if cfg.zeroPad {
			zeros = n - window.count()
		}
		y[i] = window.median(zeros)
	}
	return
}

// MedFiltOption configures MedFilt: MedFiltTruncate or MedFiltZeroPad.
type MedFiltOption func(*medFiltConfig)

// MedFiltTruncate takes the median of only the samples inside x,
// so a window near the ends is shorter, and may hold an even count.
// It's the default, as the ends stay where the signal is, wherever that is.
//
// Equivalent to:
//   - SciPy has none
//   - MATLAB medfilt1's 'truncate'
func MedFiltTruncate() MedFiltOption {
	return func(c *medFiltConfig) { c.zeroPad = false }
}

// MedFiltZeroPad fills the window past the ends with zeros,
// which pulls the ends towards 0, so it suits a signal centred on 0.
//
// Equivalent to:
//   - scipy.signal.medfilt, which has no other
//   - MATLAB medfilt1's 'zeropad', its default
func MedFiltZeroPad() MedFiltOption {
	return func(c *medFiltConfig) { c.zeroPad = true }
}

// medFiltConfig is MedFilt's settings, built from MedFiltOptions.
type medFiltConfig struct {
	zeroPad bool // pad past the ends with zeros, rather than truncate
}

// sortedWindow is a sliding window's samples, kept sorted,
// with its NaNs counted apart, as they don't sort.
type sortedWindow struct {
	sorted []float64
	nans   int
}

// add puts v into the window, in order.
func (w *sortedWindow) add(v float64) {
	if math.IsNaN(v) {
		w.nans++
		return
	}
	i, _ := slices.BinarySearch(w.sorted, v)
	w.sorted = slices.Insert(w.sorted, i, v)
}

// remove takes one sample equal to v out of the window.
func (w *sortedWindow) remove(v float64) {
	if math.IsNaN(v) {
		w.nans--
		return
	}
	i, _ := slices.BinarySearch(w.sorted, v)
	w.sorted = slices.Delete(w.sorted, i, i+1)
}

// count is how many samples the window holds, NaNs included.
func (w *sortedWindow) count() int {
	return len(w.sorted) + w.nans
}

// median is the median of the window with zeros extra zeros in it,
// as MedFiltZeroPad's padding is counted rather than stored.
// It's NaN if the window holds a NaN.
func (w *sortedWindow) median(zeros int) float64 {
	if w.nans > 0 {
		return math.NaN()
	}
	count := len(w.sorted) + zeros
	middle := count / 2
	if count%2 == 1 {
		return w.at(middle, zeros)
	}
	return floats.Midpoint(w.at(middle-1, zeros), w.at(middle, zeros))
}

// at is the k-th smallest of the window's values and zeros extra zeros.
func (w *sortedWindow) at(k, zeros int) float64 {
	negatives, _ := slices.BinarySearch(w.sorted, 0)
	switch {
	case k < negatives:
		return w.sorted[k]
	case k < negatives+zeros:
		return 0
	default:
		return w.sorted[k-zeros]
	}
}
