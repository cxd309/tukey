package conformance

import (
	"github.com/cxd309/tukey/dsp"
	"github.com/cxd309/tukey/internal/reference"
)

type filterVector struct {
	reference.Meta
	Params struct {
		B []float64 `json:"b"`
		A []float64 `json:"a"`
		X []float64 `json:"x"`
	} `json:"params"`
	Output struct {
		Y []float64 `json:"y"`
	} `json:"output"`
}

// filterSciPySuite checks dsp.Filter against scipy.signal.lfilter
// measured: FIR and simple-coefficient cases are bit-identical; Butterworth cases
// agree within ~3 ULPs (worst 6.6e-16), differing only in the order of additions;
// set ~10x above the worst for headroom across platforms
var filterSciPySuite = newSuiteOf("Filter", reference.SourceSciPy, "filter",
	reference.Tolerance{Rel: 1e-14, Scaled: 1e-14},
	filterOutputs, isWellConditioned)

// filterMATLABSuite checks dsp.Filter against MATLAB's filter
// measured: worst 5.3e-16 (2 ULP), as against SciPy; set to match
var filterMATLABSuite = newSuiteOf("Filter", reference.SourceMATLAB, "filter",
	reference.Tolerance{Rel: 1e-14, Scaled: 1e-14},
	filterOutputs, isWellConditioned)

// filterOutputs pairs dsp.Filter's output with a vector's
func filterOutputs(v filterVector) (outputs []Output, err error) {
	y, err := dsp.Filter(v.Params.B, v.Params.A, v.Params.X)
	if err != nil {
		return nil, err
	}
	outputs = []Output{{Name: "y", Got: y, Want: v.Output.Y, Scale: maxAbs(v.Params.X)}}
	return
}
