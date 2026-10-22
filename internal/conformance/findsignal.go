package conformance

import (
	"github.com/cxd309/tukey/dsp"
	"github.com/cxd309/tukey/internal/reference"
)

type findSignalVector struct {
	reference.Meta
	Params struct {
		Data   []float64 `json:"data"`
		Signal []float64 `json:"signal"`
	} `json:"params"`
	Output struct {
		Start   int       `json:"start"`
		Stop    int       `json:"stop"`
		Dist    []float64 `json:"dist"`
		Profile []float64 `json:"profile"` // SciPy only: every window's distance, for checking a whole search
	} `json:"output"`
}

// findsignal measured, on both architectures, set ~10x above the worst:
//   - direct and hybrid: exact against SciPy, and 5.9e-16 (3 ULP) against MATLAB, on delayed_walk
//   - fast: 1.5e-15 (7 ULP) against both, on one_sample on arm64
//
// Direct and hybrid against SciPy are measured relative to the distance itself, as both
// sum the closest window exactly, so a perfect match must be 0. The rest are relative to
// the window's and signal's sums of squares, as the expanded form Σd² − 2Σdt + Σt² rounds
// at that scale: MATLAB's perfect matches are up to 2^-29 rather than 0.
var (
	findSignalDirectSciPySuite  = newSuite("FindSignal (direct)", reference.SourceSciPy, "findsignal", tolerance(1e-14), findSignalOutputs(dsp.SearchDirect, false))
	findSignalDirectMATLABSuite = newSuite("FindSignal (direct)", reference.SourceMATLAB, "findsignal", tolerance(1e-14), findSignalOutputs(dsp.SearchDirect, true))
	findSignalFastSciPySuite    = newSuite("FindSignal (fast)", reference.SourceSciPy, "findsignal", tolerance(1e-13), findSignalOutputs(dsp.SearchFast, true))
	findSignalFastMATLABSuite   = newSuite("FindSignal (fast)", reference.SourceMATLAB, "findsignal", tolerance(1e-13), findSignalOutputs(dsp.SearchFast, true))
	findSignalHybridSciPySuite  = newSuite("FindSignal (hybrid)", reference.SourceSciPy, "findsignal", tolerance(1e-14), findSignalOutputs(dsp.SearchHybrid, false))
	findSignalHybridMATLABSuite = newSuite("FindSignal (hybrid)", reference.SourceMATLAB, "findsignal", tolerance(1e-14), findSignalOutputs(dsp.SearchHybrid, true))
)

// findSignalOutputs pairs dsp.FindSignal's match, by method, with a vector's.
// scaled measures the distance relative to the sums of squares rather than itself.
func findSignalOutputs(method dsp.SearchMethod, scaled bool) func(v findSignalVector) ([]Output, error) {
	return func(v findSignalVector) (outputs []Output, err error) {
		match, err := dsp.FindSignal(v.Params.Data, v.Params.Signal, dsp.SignalSearch(method))
		if err != nil {
			return nil, err
		}
		var scale float64
		if scaled {
			scale = sumSquares(v.Params.Signal) + sumSquares(v.Params.Data[v.Output.Start:v.Output.Stop+1])
		}
		outputs = []Output{
			{Name: "start", Got: []float64{float64(match.Start)}, Want: []float64{float64(v.Output.Start)}},
			{Name: "stop", Got: []float64{float64(match.Stop)}, Want: []float64{float64(v.Output.Stop)}},
			{Name: "dist", Got: []float64{match.Distance}, Want: v.Output.Dist, Scale: scale},
		}
		return
	}
}

// sumSquares is the sum of the squares of x.
func sumSquares(x []float64) (sum float64) {
	for _, v := range x {
		sum += v * v
	}
	return
}
