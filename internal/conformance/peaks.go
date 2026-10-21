package conformance

import (
	"github.com/cxd309/tukey/dsp"
	"github.com/cxd309/tukey/internal/reference"
)

type findPeaksVector struct {
	reference.Meta
	Params struct {
		X          []float64 `json:"x"`
		Height     *float64  `json:"height"`
		Threshold  *float64  `json:"threshold"`
		Distance   *int      `json:"distance"`
		Prominence *float64  `json:"prominence"`
		Width      *float64  `json:"width"`
	} `json:"params"`
	Output struct {
		Peaks   []int     `json:"peaks"`
		Heights []float64 `json:"heights"`
	} `json:"output"`
}

// find_peaks is exact against both. The peaks are indices and the heights are
// copied from x, so there's no rounding to allow for.
var (
	findPeaksSciPySuite  = newSuite("FindPeaks", reference.SourceSciPy, "find_peaks", reference.Tolerance{}, findPeaksSciPyOutputs)
	findPeaksMATLABSuite = newSuite("FindPeaks", reference.SourceMATLAB, "find_peaks", reference.Tolerance{}, findPeaksMATLABOutputs)
)

// findPeaksSciPyOutputs runs dsp.FindPeaks with the vector's options and
// tukey's defaults, which are SciPy's rules.
func findPeaksSciPyOutputs(v findPeaksVector) (outputs []Output, err error) {
	outputs, err = findPeaksOutputs(v, peakOptions(v, 0))
	return
}

// findPeaksMATLABOutputs runs dsp.FindPeaks with findpeaks' rules. MATLAB's
// vectors keep SciPy's params, and findpeaks' MinPeakDistance d keeps peaks
// more than d apart, which is MinPeakDistance(d+1).
func findPeaksMATLABOutputs(v findPeaksVector) (outputs []Output, err error) {
	opts := append(peakOptions(v, 1), dsp.MATLABPeaks())
	outputs, err = findPeaksOutputs(v, opts)
	return
}

// peakOptions is the vector's params as PeakOptions, with distanceOffset
// added to its distance.
func peakOptions(v findPeaksVector, distanceOffset int) (opts []dsp.PeakOption) {
	p := v.Params
	if p.Height != nil {
		opts = append(opts, dsp.MinPeakHeight(*p.Height))
	}
	if p.Threshold != nil {
		opts = append(opts, dsp.PeakThreshold(*p.Threshold))
	}
	if p.Distance != nil {
		opts = append(opts, dsp.MinPeakDistance(*p.Distance+distanceOffset))
	}
	if p.Prominence != nil {
		opts = append(opts, dsp.MinPeakProminence(*p.Prominence))
	}
	if p.Width != nil {
		opts = append(opts, dsp.MinPeakWidth(*p.Width))
	}
	return
}

// findPeaksOutputs pairs dsp.FindPeaks' peaks and heights with a vector's.
func findPeaksOutputs(v findPeaksVector, opts []dsp.PeakOption) (outputs []Output, err error) {
	peaks, err := dsp.FindPeaks(v.Params.X, opts...)
	if err != nil {
		return nil, err
	}
	indices := make([]float64, len(peaks))
	heights := make([]float64, len(peaks))
	for k, p := range peaks {
		indices[k], heights[k] = float64(p.Index), p.Height
	}
	outputs = []Output{
		{Name: "peaks", Got: indices, Want: indexFloats(v.Output.Peaks)},
		{Name: "heights", Got: heights, Want: v.Output.Heights},
	}
	return
}

type peakProminencesVector struct {
	reference.Meta
	Params struct {
		X     []float64 `json:"x"`
		Peaks []int     `json:"peaks"`
	} `json:"params"`
	Output struct {
		Prominences []float64 `json:"prominences"`
		LeftBases   []int     `json:"left_bases"`  // SciPy only
		RightBases  []int     `json:"right_bases"` // SciPy only
	} `json:"output"`
}

// peak_prominences measured: exact against both, on both architectures.
// Each prominence is one subtraction of two samples of x, the same in all three,
// so exact is expected, but the tolerance stays at 1e-14.
var (
	peakProminencesSciPySuite  = newSuite("PeakProminences", reference.SourceSciPy, "peak_prominences", tolerance(1e-14), peakProminencesOutputs)
	peakProminencesMATLABSuite = newSuite("PeakProminences", reference.SourceMATLAB, "peak_prominences", tolerance(1e-14), peakProminencesOutputs)
)

// peakProminencesOutputs pairs dsp.PeakProminences' output with a vector's,
// and its bases too where the vector has them.
func peakProminencesOutputs(v peakProminencesVector) (outputs []Output, err error) {
	got, err := dsp.PeakProminences(v.Params.X, v.Params.Peaks)
	if err != nil {
		return nil, err
	}
	prominences := make([]float64, len(got))
	left := make([]float64, len(got))
	right := make([]float64, len(got))
	for k, p := range got {
		prominences[k], left[k], right[k] = p.Prominence, float64(p.LeftBase), float64(p.RightBase)
	}
	outputs = []Output{{Name: "prominences", Got: prominences, Want: v.Output.Prominences}}
	if v.Output.LeftBases != nil {
		outputs = append(outputs,
			Output{Name: "left_bases", Got: left, Want: indexFloats(v.Output.LeftBases)},
			Output{Name: "right_bases", Got: right, Want: indexFloats(v.Output.RightBases)},
		)
	}
	return
}

type peakWidthsVector struct {
	reference.Meta
	Params struct {
		X         []float64 `json:"x"`
		Peaks     []int     `json:"peaks"`
		RelHeight float64   `json:"rel_height"`
	} `json:"params"`
	Output struct {
		Widths       []float64 `json:"widths"`
		WidthHeights []float64 `json:"width_heights"` // SciPy only
		LeftIPs      []float64 `json:"left_ips"`      // SciPy only
		RightIPs     []float64 `json:"right_ips"`     // SciPy only
	} `json:"output"`
}

// peak_widths measured: exact against SciPy, and worst 2.7e-15 (12 ULP) against
// MATLAB, on abs_noise, on both architectures. Set ~10x above.
var (
	peakWidthsSciPySuite  = newSuite("PeakWidths", reference.SourceSciPy, "peak_widths", tolerance(1e-14), peakWidthsOutputs)
	peakWidthsMATLABSuite = newSuite("PeakWidths", reference.SourceMATLAB, "peak_widths", tolerance(1e-13), peakWidthsOutputs)
)

// peakWidthsOutputs pairs dsp.PeakWidths' output with a vector's, and the
// heights and crossings too where the vector has them.
func peakWidthsOutputs(v peakWidthsVector) (outputs []Output, err error) {
	got, err := dsp.PeakWidths(v.Params.X, v.Params.Peaks, v.Params.RelHeight)
	if err != nil {
		return nil, err
	}
	widths := make([]float64, len(got))
	heights := make([]float64, len(got))
	left := make([]float64, len(got))
	right := make([]float64, len(got))
	for k, w := range got {
		widths[k], heights[k], left[k], right[k] = w.Width, w.Height, w.Left, w.Right
	}
	outputs = []Output{{Name: "widths", Got: widths, Want: v.Output.Widths}}
	if v.Output.WidthHeights != nil {
		outputs = append(outputs,
			Output{Name: "width_heights", Got: heights, Want: v.Output.WidthHeights},
			Output{Name: "left_ips", Got: left, Want: v.Output.LeftIPs},
			Output{Name: "right_ips", Got: right, Want: v.Output.RightIPs},
		)
	}
	return
}

// indexFloats is indices as float64, to compare as an Output.
func indexFloats(indices []int) (floats []float64) {
	floats = make([]float64, len(indices))
	for k, i := range indices {
		floats[k] = float64(i)
	}
	return
}
