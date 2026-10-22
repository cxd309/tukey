package dsp

import (
	"cmp"
	"fmt"
	"math"
	"slices"
)

// Peak is a local maximum of a signal: its index, and the signal's value there
type Peak struct {
	Index  int
	Height float64
}

// PeakProminence is how far a peak stands above the higher of the lowest points
// either side of it, before the signal rises above the peak again.
// LeftBase and RightBase are the indices of those lowest points.
type PeakProminence struct {
	Prominence float64
	LeftBase   int
	RightBase  int
}

// PeakWidth is a peak's width at Height,
// measured between where the signal crosses Height either side of the peak.
// Left and Right are those crossings, as fractional indicies
// interpolated between samples.
type PeakWidth struct {
	Width  float64
	Height float64
	Left   float64
	Right  float64
}

// PlateauPosition is the sample FindPeaks reports for a flat-topped peak
type PlateauPosition int

const (
	// PlateauMiddle reports the middle sample,
	// rounding down for an even plateau, default same as SciPy
	PlateauMiddle PlateauPosition = iota
	// PlateauFirst reports the first sample, same as MATLAB
	PlateauFirst
)

// TieBreak is which of two equal peaks MinDistance keeps
type TieBreak int

const (
	// TieLater keeps the later peak, default same as SciPy
	TieLater TieBreak = iota
	// TieEarlier keeps the earlier peak, same as MATLAB
	TieEarlier
)

// PeakOption configures FindPeaks.
// Options are applied in order, so a later option overrides an earlier one
type PeakOption func(*peakConfig)

// peakConfig is FindPeak' settings, built from PeakOptions.
// A limit of -Inf, or a distance of 1, means no limit
type peakConfig struct {
	minHeight     float64
	threshold     float64
	minDistance   int
	minProminence float64
	minWidth      float64
	plateau       PlateauPosition
	ties          TieBreak
	distanceLast  bool
}

// MinPeakHeight keeps peaks at least h high
//
// Equivalent to SciPy height=h,  MATLAB "MinPeakHeight"
func MinPeakHeight(h float64) PeakOption {
	return func(c *peakConfig) { c.minHeight = h }
}

// Treshold keeps peaks at least t above both of their neighbouring samples.
//
// Equivalent to SciPy threshold=t, MATLAB "PeakThreshold"
func PeakThreshold(t float64) PeakOption {
	return func(c *peakConfig) { c.threshold = t }

}

// MinPeakDistance keeps peaks at least d samples apart,
// removing the lower of any two that are closer, highest peak first.
//
// Equivalent of SciPy distance=d.
// MATLAB "MinPeakDistance", d keeps peaks more than d apart, which is MinPeakDistance(d+1)
func MinPeakDistance(d int) PeakOption {
	return func(c *peakConfig) { c.minDistance = d }
}

// MinPeakProminence keeps peaks with a prominence of at least p (see PeakProminences)
//
// Equivalent to SciPy prominence=p, MATLAB "MinPeakProminence"
func MinPeakProminence(p float64) PeakOption {
	return func(c *peakConfig) { c.minProminence = p }
}

// MinPeakWidth keeps peaks at least w samples wide at half their prominence (see PeakWidths)
//
// Equivalent to SciPy width=w, MATLAB "MinPeakWidth"
func MinPeakWidth(w float64) PeakOption {
	return func(c *peakConfig) { c.minWidth = w }
}

// Plateau sets which sample of a flat-topped peak is reported
func Plateau(p PlateauPosition) PeakOption {
	return func(c *peakConfig) { c.plateau = p }
}

// Ties sets which of two equal peaks MinDistance keeps
func Ties(t TieBreak) PeakOption {
	return func(c *peakConfig) { c.ties = t }
}

// DistanceLast applies MinDistance after MinProminence and MinWidth,
// so a peak that fails those can't remove a neighbour that passes them.
// By default distance is applied first, as SciPy.
func DistanceLast() PeakOption {
	return func(c *peakConfig) { c.distanceLast = true }
}

// MATLABPeaks applies MATLAB findpeaks' rules:
//   - first sample of a plateau
//   - ties to the earlier peak
//   - distance applied last
//
// It doesn't change MinDistance, which is d+1 for MATLAB d
func MATLABPeaks() PeakOption {
	return func(c *peakConfig) {
		c.plateau = PlateauFirst
		c.ties = TieEarlier
		c.distanceLast = true
	}
}

// FindPeaks returns the local maxima of x in index order.
// A local maximum is a sample, of a flat run of samples, higher than the
// samples either side.
// The first and last samples are never peaks, so x shorter than 3 has none.
//
// Options narrow the peaks down. Each of the following just a peak on its own:
//   - MinHeight
//   - Threshold
//   - MinProminence
//   - MinWidth
//
// MinDistance removes the lower of any two peaks that are too close together.
//
// Returns ErrInvalidParameter for a MinDistance below 1.
//
// Equivalent to scipy.signal.find_peaks.
// With MATLABPeaks, equivalent to MATLAB findpeaks,
// which errors for x shorter than 3 and for a MinPeakDIstance of len(x)-1 or more
func FindPeaks(x []float64, opts ...PeakOption) (peaks []Peak, err error) {
	cfg, err := newPeakConfig(opts)
	if err != nil {
		return
	}
	indices := localMaxima(x, cfg.plateau)
	indices = slices.DeleteFunc(indices, func(i int) bool {
		return x[i] < cfg.minHeight || min(x[i]-x[i-1], x[i]-x[i+1]) < cfg.threshold
	})
	if !cfg.distanceLast {
		indices = byDistance(x, indices, cfg.minDistance, cfg.ties)
	}
	indices = byShape(x, indices, cfg)
	if cfg.distanceLast {
		indices = byDistance(x, indices, cfg.minDistance, cfg.ties)
	}
	peaks = make([]Peak, len(indices))
	for k, i := range indices {
		peaks[k] = Peak{Index: i, Height: x[i]}
	}
	return
}

// PeakProminences returns the prominence of each of peaks in x,
// and the bases it's measured from.
//
// Returns ErrInvalidParameter for a peak outside x.
//
// Equivalent to scipy.signal.peak_prominences
// Equivalent to the fourth output of MATLAB findpeaks
func PeakProminences(x []float64, peaks []int) (prominences []PeakProminence, err error) {
	if err = checkPeaks(x, peaks); err != nil {
		return
	}
	prominences = make([]PeakProminence, len(peaks))
	for k, p := range peaks {
		prominences[k] = peakProminence(x, p)
	}
	return
}

// PeakWidths returns the width of each of peaks in x,
// measured relHeight of the way down its prominence:
// 0.5 is half the prominence, 1 is at its bases.
//
// Returns ErrInvalidParameter for a negative relHeight
// or a peak outside x
//
// Equivalent to scipy.signal.peak_widths.
// Equivalent to third output of MATLAB findpeaks at relHeight 0.5
func PeakWidths(x []float64, peaks []int, relHeight float64) (widths []PeakWidth, err error) {
	if !(relHeight >= 0) {
		return nil, fmt.Errorf("%w: relHeight %g, must be at least 0", ErrInvalidParameter, relHeight)
	}
	if err = checkPeaks(x, peaks); err != nil {
		return
	}
	widths = make([]PeakWidth, len(peaks))
	for k, p := range peaks {
		widths[k] = peakWidth(x, p, relHeight, peakProminence(x, p))
	}
	return
}

// newPeakConfig builds FindPeaks' settings from opts
func newPeakConfig(opts []PeakOption) (cfg peakConfig, err error) {
	cfg = peakConfig{
		minHeight:     math.Inf(-1),
		threshold:     math.Inf(-1),
		minDistance:   1,
		minProminence: math.Inf(-1),
		minWidth:      math.Inf(-1),
	}
	for _, opt := range opts {
		opt(&cfg)
	}
	if cfg.minDistance < 1 {
		err = fmt.Errorf("%w: MinDistance %d, must be at least 1", ErrInvalidParameter, cfg.minDistance)
	}
	return
}

// localMaxima is the index of every local maximum of x, with a plateau
// reported at the sample plateau says
func localMaxima(x []float64, plateau PlateauPosition) (indices []int) {
	n := len(x)
	for i := 1; i < n-1; i++ {
		if !(x[i-1] < x[i]) {
			continue
		}
		ahead := i + 1 // the first sample after the plateau starting at i
		for ahead < n-1 && x[ahead] == x[i] {
			ahead++
		}
		if x[ahead] < x[i] {
			indices = append(indices, plateau.pick(i, ahead-1))
			i = ahead
		}
	}
	return
}

// pick is the sample reported for a plateau from first to last
func (p PlateauPosition) pick(first, last int) (index int) {
	if p == PlateauFirst {
		return first
	}
	return (first + last) / 2
}

// byDistance removes the lower of any two of indices less than distance
// apart, working down from the highest peak, with ties broken by ties
func byDistance(x []float64, indices []int, distance int, ties TieBreak) (kept []int) {
	if distance <= 1 {
		return indices
	}
	order := make([]int, len(indices)) // positions in indices, highest peak first
	for k := range order {
		order[k] = k
	}
	slices.SortFunc(order, func(a, b int) int {
		if c := cmp.Compare(x[indices[b]], x[indices[a]]); c != 0 {
			return c
		}
		if ties == TieEarlier {
			return cmp.Compare(a, b)
		}
		return cmp.Compare(b, a)
	})
	removed := make([]bool, len(indices))
	for _, j := range order {
		if removed[j] {
			continue
		}
		for k := j - 1; k >= 0 && indices[j]-indices[k] < distance; k-- {
			removed[k] = true
		}
		for k := j + 1; k < len(indices) && indices[k]-indices[j] < distance; k++ {
			removed[k] = true
		}
	}
	for k, i := range indices {
		if !removed[k] {
			kept = append(kept, i)
		}
	}
	return
}

// byShape removes the indices that fail cfg's MinProminence or MinWidth.
// Each peak is judged on the whole of x, so the order doesn't matter here
func byShape(x []float64, indices []int, cfg peakConfig) (kept []int) {
	byProminence := cfg.minProminence > math.Inf(-1)
	byWidth := cfg.minWidth > math.Inf(-1)
	if !byProminence && !byWidth {
		return indices
	}
	kept = slices.DeleteFunc(indices, func(i int) bool {
		prominence := peakProminence(x, i)
		if prominence.Prominence < cfg.minProminence {
			return true
		}
		return byWidth && peakWidth(x, i, 0.5, prominence).Width < cfg.minWidth
	})
	return
}

// peakProminence walks out from peak each way until x rises above it,
// and measures from the higher of the lowest points found
func peakProminence(x []float64, peak int) (p PeakProminence) {
	height := x[peak]
	leftMin, rightMin := height, height
	p.LeftBase, p.RightBase = peak, peak
	for i := peak; i >= 0 && x[i] <= height; i-- {
		if x[i] < leftMin {
			leftMin, p.LeftBase = x[i], i
		}
	}
	for i := peak; i < len(x) && x[i] <= height; i++ {
		if x[i] < rightMin {
			rightMin, p.RightBase = x[i], i
		}
	}
	p.Prominence = height - max(leftMin, rightMin)
	return
}

// peakWidth walks out from peak each way, no further than its bases,
// until x drops to the width height, then interpolates the crossing between samples
func peakWidth(x []float64, peak int, relHeight float64, prominence PeakProminence) (w PeakWidth) {
	w.Height = x[peak] - prominence.Prominence*relHeight
	i := peak
	for prominence.LeftBase < i && w.Height < x[i] {
		i--
	}
	w.Left = float64(i)
	if x[i] < w.Height {
		w.Left += (w.Height - x[i]) / (x[i+1] - x[i])
	}
	i = peak
	for i < prominence.RightBase && w.Height < x[i] {
		i++
	}
	w.Right = float64(i)
	if x[i] < w.Height {
		w.Right -= (w.Height - x[i]) / (x[i-1] - x[i])
	}
	w.Width = w.Right - w.Left
	return
}

// checkPeaks returns ErrInvalidParameter for any of peaks outside x
func checkPeaks(x []float64, peaks []int) (err error) {
	for _, p := range peaks {
		if p < 0 || p >= len(x) {
			return fmt.Errorf("%w: peak %d outside x of length %d", ErrInvalidParameter, p, len(x))
		}
	}
	return
}
