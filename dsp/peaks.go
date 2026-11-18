package dsp

import (
	"cmp"
	"fmt"
	"math"
	"slices"
)

// FindPeaks returns the local maxima of x, in index order.
// A local maximum is a sample, or a flat run of samples, higher than the samples either side.
// The first and last samples are never peaks, so an x shorter than 3 has none.
// Options narrow the peaks down:
//   - MinPeakHeight, PeakThreshold, MinPeakProminence and MinPeakWidth each judge a peak on its own
//   - MinPeakDistance removes the lower of any two peaks that are too close together
//
// A NaN is never a peak, nor is a sample beside one, as it isn't known to be higher.
// A peak whose prominence or width search reaches a NaN has neither, so fails MinPeakProminence and MinPeakWidth.
// An infinite limit is allowed: -Inf is no limit, and +Inf passes no peak.
//
//	peaks, err := dsp.FindPeaks(x, dsp.MinPeakProminence(0.5), dsp.MinPeakDistance(10))
//
// Units: Peak's Index is a sample index, and its Height is in x's unit
//
// Returns: ErrInvalidParameter for a MinPeakDistance below 1, or a NaN limit.
//
// Equivalent to:
//   - scipy.signal.find_peaks(x)
//   - MATLAB findpeaks(x), with MATLABPeaks, which errors for an x shorter than 3,
//     and for a 'MinPeakDistance' of len(x)-1 or more
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

// PeakProminences returns the prominence of each of peaks in x, and the bases it's measured from:
// how far the peak stands above the higher of the lowest points either side of it,
// before the signal rises above the peak again.
// A peak whose search reaches a NaN before a higher sample has a NaN prominence, as the NaN might be higher.
//
//	prominences, err := dsp.PeakProminences(x, []int{12, 40})
//
// Units:
//   - peaks are sample indices
//   - each Prominence is in x's unit, and LeftBase and RightBase are sample indices
//
// Returns: ErrInvalidParameter for a peak outside x.
//
// Equivalent to:
//   - scipy.signal.peak_prominences(x, peaks)
//   - MATLAB findpeaks' fourth output
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

// PeakWidths returns the width of each of peaks in x, measured relHeight of the way down its prominence:
// 0.5 is half the prominence, and 1 is at its bases.
// The crossings are interpolated between samples.
// A peak with a NaN prominence, from PeakProminences, has a NaN width.
//
//	widths, err := dsp.PeakWidths(x, []int{12, 40}, 0.5)
//
// Units:
//   - peaks are sample indices
//   - relHeight is a fraction of the prominence
//   - each Width is in samples, its Height in x's unit, and Left and Right are fractional sample indices
//
// Returns: ErrInvalidParameter for a relHeight that isn't at least 0 and finite, or a peak outside x.
//
// Equivalent to:
//   - scipy.signal.peak_widths(x, peaks, rel_height=relHeight)
//   - MATLAB findpeaks' third output, at relHeight 0.5
func PeakWidths(x []float64, peaks []int, relHeight float64) (widths []PeakWidth, err error) {
	if !(relHeight >= 0) || math.IsInf(relHeight, 1) {
		return nil, fmt.Errorf("%w: relHeight %g, must be at least 0 and finite", ErrInvalidParameter, relHeight)
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

// Peak is a local maximum of a signal: its index, and the signal's value there.
type Peak struct {
	Index  int     // a sample index
	Height float64 // in the signal's unit
}

// PeakProminence is how far a peak stands above the higher of the lowest points either side of it,
// before the signal rises above the peak again.
// LeftBase and RightBase are the indices of those lowest points.
type PeakProminence struct {
	Prominence float64 // in the signal's unit
	LeftBase   int     // a sample index
	RightBase  int     // a sample index
}

// PeakWidth is a peak's width at Height,
// measured between where the signal crosses Height either side of the peak.
// Left and Right are those crossings, as fractional indices interpolated between samples.
type PeakWidth struct {
	Width  float64 // in samples, Right - Left
	Height float64 // in the signal's unit
	Left   float64 // a fractional sample index
	Right  float64 // a fractional sample index
}

// PeakOption configures FindPeaks.
// Options are applied in order, so a later one overrides an earlier one.
type PeakOption func(*peakConfig)

// MinPeakHeight keeps peaks at least h high.
// -Inf is no limit, and +Inf passes no peak.
//
// Units: h is in x's unit
//
// Equivalent to:
//   - scipy.signal.find_peaks' height=h
//   - MATLAB findpeaks' 'MinPeakHeight', h
func MinPeakHeight(h float64) PeakOption {
	return func(c *peakConfig) { c.minHeight = h }
}

// PeakThreshold keeps peaks at least t above both of their neighbouring samples.
// -Inf is no limit, and +Inf passes no peak.
//
// Units: t is in x's unit
//
// Equivalent to:
//   - scipy.signal.find_peaks' threshold=t
//   - MATLAB findpeaks' 'Threshold', t
func PeakThreshold(t float64) PeakOption {
	return func(c *peakConfig) { c.threshold = t }

}

// MinPeakDistance keeps peaks at least d samples apart,
// removing the lower of any two that are closer, highest peak first.
//
// Units: d is a number of samples
//
// Equivalent to:
//   - scipy.signal.find_peaks' distance=d
//   - MATLAB findpeaks' 'MinPeakDistance', d-1, as MATLAB keeps peaks more than its distance apart
func MinPeakDistance(d int) PeakOption {
	return func(c *peakConfig) { c.minDistance = d }
}

// MinPeakProminence keeps peaks with a prominence, from PeakProminences, of at least p.
// -Inf is no limit, and +Inf passes no peak.
//
// Units: p is in x's unit
//
// Equivalent to:
//   - scipy.signal.find_peaks' prominence=p
//   - MATLAB findpeaks' 'MinPeakProminence', p
func MinPeakProminence(p float64) PeakOption {
	return func(c *peakConfig) { c.minProminence = p }
}

// MinPeakWidth keeps peaks at least w samples wide at half their prominence, from PeakWidths.
// -Inf is no limit, and +Inf passes no peak.
//
// Units: w is in samples
//
// Equivalent to:
//   - scipy.signal.find_peaks' width=w
//   - MATLAB findpeaks' 'MinPeakWidth', w
func MinPeakWidth(w float64) PeakOption {
	return func(c *peakConfig) { c.minWidth = w }
}

// PlateauMiddle reports a flat-topped peak at its middle sample, rounding down for an even plateau.
// It's the default.
//
// Equivalent to:
//   - scipy.signal.find_peaks, which has no other
//   - MATLAB has none
func PlateauMiddle() PeakOption {
	return func(c *peakConfig) { c.plateau = plateauMiddle }
}

// PlateauFirst reports a flat-topped peak at its first sample.
//
// Equivalent to:
//   - SciPy has none
//   - MATLAB findpeaks, which has no other
func PlateauFirst() PeakOption {
	return func(c *peakConfig) { c.plateau = plateauFirst }
}

// TieLater keeps the later of two equal peaks closer than MinPeakDistance.
// It's the default.
//
// Equivalent to:
//   - scipy.signal.find_peaks, which has no other
//   - MATLAB has none
func TieLater() PeakOption {
	return func(c *peakConfig) { c.ties = tieLater }
}

// TieEarlier keeps the earlier of two equal peaks closer than MinPeakDistance.
//
// Equivalent to:
//   - SciPy has none
//   - MATLAB findpeaks, which has no other
func TieEarlier() PeakOption {
	return func(c *peakConfig) { c.ties = tieEarlier }
}

// DistanceLast applies MinPeakDistance after MinPeakProminence and MinPeakWidth,
// so a peak that fails those can't remove a neighbour that passes them.
// Without it, distance is applied first.
//
// Equivalent to:
//   - SciPy has none, as find_peaks applies distance first
//   - MATLAB findpeaks, which has no other
func DistanceLast() PeakOption {
	return func(c *peakConfig) { c.distanceLast = true }
}

// MATLABPeaks applies MATLAB findpeaks' rules, same as PlateauFirst, TieEarlier and DistanceLast together:
//   - a plateau's first sample
//   - ties to the earlier peak
//   - distance applied last
//
// It doesn't change MinPeakDistance, which is d+1 for MATLAB's d.
//
// Equivalent to:
//   - SciPy has none
//   - MATLAB findpeaks
func MATLABPeaks() PeakOption {
	return func(c *peakConfig) {
		c.plateau = plateauFirst
		c.ties = tieEarlier
		c.distanceLast = true
	}
}

// peakConfig is FindPeaks' settings, built from PeakOptions.
// A limit of -Inf, or a distance of 1, means no limit.
type peakConfig struct {
	minHeight     float64
	threshold     float64
	minDistance   int
	minProminence float64
	minWidth      float64
	plateau       plateauPosition
	ties          tieBreak
	distanceLast  bool
}

// plateauPosition is the sample FindPeaks reports for a flat-topped peak.
type plateauPosition int

const (
	plateauMiddle plateauPosition = iota
	plateauFirst
)

// tieBreak is which of two equal peaks MinPeakDistance keeps.
type tieBreak int

const (
	tieLater tieBreak = iota
	tieEarlier
)

// newPeakConfig builds FindPeaks' settings from opts, returning ErrInvalidParameter for a NaN limit or a distance below 1.
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
	limits := []struct {
		name  string
		value float64
	}{
		{"MinPeakHeight", cfg.minHeight},
		{"PeakThreshold", cfg.threshold},
		{"MinPeakProminence", cfg.minProminence},
		{"MinPeakWidth", cfg.minWidth},
	}
	for _, limit := range limits {
		if math.IsNaN(limit.value) {
			return cfg, fmt.Errorf("%w: %s is NaN", ErrInvalidParameter, limit.name)
		}
	}
	if cfg.minDistance < 1 {
		err = fmt.Errorf("%w: MinDistance %d, must be at least 1", ErrInvalidParameter, cfg.minDistance)
	}
	return
}

// localMaxima is the index of every local maximum of x, with a plateau reported at the sample plateau says.
func localMaxima(x []float64, plateau plateauPosition) (indices []int) {
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

// pick is the sample reported for a plateau from first to last.
func (p plateauPosition) pick(first, last int) (index int) {
	if p == plateauFirst {
		return first
	}
	return (first + last) / 2
}

// byDistance removes the lower of any two of indices less than distance apart,
// working down from the highest peak, with ties broken by ties.
func byDistance(x []float64, indices []int, distance int, ties tieBreak) (kept []int) {
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
		if ties == tieEarlier {
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

// byShape removes the indices that fail cfg's MinPeakProminence or MinPeakWidth.
// Each peak is judged on the whole of x, so the order doesn't matter here.
func byShape(x []float64, indices []int, cfg peakConfig) (kept []int) {
	byProminence := cfg.minProminence > math.Inf(-1)
	byWidth := cfg.minWidth > math.Inf(-1)
	if !byProminence && !byWidth {
		return indices
	}
	kept = slices.DeleteFunc(indices, func(i int) bool {
		prominence := peakProminence(x, i)
		if !(prominence.Prominence >= cfg.minProminence) { // NaN fails too
			return true
		}
		return byWidth && !(peakWidth(x, i, 0.5, prominence).Width >= cfg.minWidth)
	})
	return
}

// peakProminence walks out from peak each way until x rises above it,
// and measures from the higher of the lowest points found.
// A NaN met on the way leaves the prominence unknown, so NaN.
func peakProminence(x []float64, peak int) (p PeakProminence) {
	height := x[peak]
	leftMin, rightMin := height, height
	p.LeftBase, p.RightBase = peak, peak
	i := peak
	for ; i >= 0 && x[i] <= height; i-- {
		if x[i] < leftMin {
			leftMin, p.LeftBase = x[i], i
		}
	}
	metNaN := i >= 0 && math.IsNaN(x[i])
	for i = peak; i < len(x) && x[i] <= height; i++ {
		if x[i] < rightMin {
			rightMin, p.RightBase = x[i], i
		}
	}
	metNaN = metNaN || i < len(x) && math.IsNaN(x[i])
	p.Prominence = height - max(leftMin, rightMin)
	if metNaN {
		p.Prominence = math.NaN()
	}
	return
}

// peakWidth walks out from peak each way, no further than its bases, until x drops to the width's height,
// then interpolates the crossing between samples.
func peakWidth(x []float64, peak int, relHeight float64, prominence PeakProminence) (w PeakWidth) {
	w.Height = x[peak] - prominence.Prominence*relHeight
	if math.IsNaN(prominence.Prominence) {
		w.Width, w.Left, w.Right = math.NaN(), math.NaN(), math.NaN()
		return
	}
	i := peak
	for prominence.LeftBase < i && w.Height < x[i] {
		i--
	}
	w.Left = float64(i)
	switch {
	case math.IsInf(x[i], -1): // the crossing's limit is the next sample
		w.Left++
	case x[i] < w.Height:
		w.Left += (w.Height - x[i]) / (x[i+1] - x[i])
	}
	i = peak
	for i < prominence.RightBase && w.Height < x[i] {
		i++
	}
	w.Right = float64(i)
	switch {
	case math.IsInf(x[i], -1):
		w.Right--
	case x[i] < w.Height:
		w.Right -= (w.Height - x[i]) / (x[i-1] - x[i])
	}
	w.Width = w.Right - w.Left
	return
}

// checkPeaks returns ErrInvalidParameter for any of peaks outside x.
func checkPeaks(x []float64, peaks []int) (err error) {
	for _, p := range peaks {
		if p < 0 || p >= len(x) {
			return fmt.Errorf("%w: peak %d outside x of length %d", ErrInvalidParameter, p, len(x))
		}
	}
	return
}
