package dsp

import (
	"fmt"
	"math"
	"slices"
)

// MedFiltEdge is how MedFilt treats a window that runs past the ends of x
type MedFiltEdge int

const (
	// MedFiltTruncate takes the median of only the samples inside x,
	// so a window near the ends is shorter, and uneven.
	// The default: the ends stay where the signal is, wherever that is
	MedFiltTruncate MedFiltEdge = iota
	// MedFiltZeroPad fills the window past the ends with zeros,
	// which pulls the ends towards 0, it suits a signal centred on 0
	MedFiltZeroPad
)

// MedFiltOption configures MedFilt
type MedFiltOption func(*medFiltConfig)

// medFiltConfig is MedFilt's settings, built from MedFiltOptions
type medFiltConfig struct {
	edge MedFiltEdge
}

// MedFiltEdges sets how MedFilt treats the ends. MedFiltTruncate by default.
func MedFiltEdges(edge MedFiltEdge) MedFiltOption {
	return func(c *medFiltConfig) { c.edge = edge }
}

// MedFilt replaces each sample of x with the median of the n samples centred on it,
// which removes spikes without blurring steps. n must be odd.
// A window holding a NaN gives NaN.
// A window with an even count, which only MedFiltTruncate gives near the ends,
// takes the midpoint of its two middle values, rounded once, so y is exact.
//
// It keeps the window sorted as it slides, O(n) per sample.
//
// Returns:
//   - ErrInvalidLength for an n that isn't odd and positive
//   - ErrInvalidMode for an unknown MedFiltEdge
//
// Equivalent to:
//   - matlab medfilt1(x, n, 'truncate'), or medfilt1(x, n) with MedFiltZeroPad
//   - scipy.signal.medfilt(x, n) with MedFiltZeroPad, except at NaN, where
//     its result depends on how NaN sorts
func MedFilt(x []float64, n int, opts ...MedFiltOption) (y []float64, err error) {
	var cfg medFiltConfig
	for _, opt := range opts {
		opt(&cfg)
	}
	switch {
	case n < 1 || n%2 == 0:
		return nil, fmt.Errorf("%w: MedFilt needs an odd, positive window length, got %d", ErrInvalidLength, n)
	case cfg.edge != MedFiltTruncate && cfg.edge != MedFiltZeroPad:
		return nil, fmt.Errorf("%w: MedFiltEdge %d", ErrInvalidMode, cfg.edge)
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
		if cfg.edge == MedFiltZeroPad {
			zeros = n - window.count()
		}
		y[i] = window.median(zeros)
	}
	return
}

// sortedWindow is a sliding window's samples, kept sorted,
// with its NaNs counted apart, as they don't sort
type sortedWindow struct {
	sorted []float64
	nans   int
}

// count is how many samples the window holds, NaN included
func (w *sortedWindow) count() int {
	return len(w.sorted) + w.nans
}

// add puts v into the window, in order
func (w *sortedWindow) add(v float64) {
	if math.IsNaN(v) {
		w.nans++
		return
	}
	i, _ := slices.BinarySearch(w.sorted, v)
	w.sorted = slices.Insert(w.sorted, i, v)
}

// remove takes one sample equal to v out of the window
func (w *sortedWindow) remove(v float64) {
	if math.IsNaN(v) {
		w.nans--
		return
	}
	i, _ := slices.BinarySearch(w.sorted, v)
	w.sorted = slices.Delete(w.sorted, i, i+1)
}

// median is the median of the window with zeros extra zeros in it,
// as MedFiltZeroPad's padding is counted rather than stored:
// NaN if the window holds a NaN
func (w *sortedWindow) median(zeros int) float64 {
	if w.nans > 0 {
		return math.NaN()
	}
	count := len(w.sorted) + zeros
	middle := count / 2
	if count%2 == 1 {
		return w.at(middle, zeros)
	}
	return midpoint(w.at(middle-1, zeros), w.at(middle, zeros))
}

// at is the k-th smallest of the window's values and zeros extra zeros
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

// midpoint is (a + b)/2 rounded once.
// a + b rounds once, and halving it is exact unless the half is subnormal,
// and then a + b, on the subnormals' grid, was exact.
// Only an overflowing a + b needs a/2 + b/2, where halving is exact.
func midpoint(a, b float64) float64 {
	sum := a + b
	if math.IsInf(sum, 0) && !math.IsInf(a, 0) && !math.IsInf(b, 0) {
		return a/2 + b/2
	}
	return sum / 2
}
