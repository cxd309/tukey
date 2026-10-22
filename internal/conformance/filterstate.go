package conformance

import (
	"github.com/cxd309/tukey/dsp"
	"github.com/cxd309/tukey/internal/reference"
)

type lfilterStateVector struct {
	reference.Meta
	Params struct {
		B  []float64 `json:"b"`
		A  []float64 `json:"a"`
		X  []float64 `json:"x"`
		Zi []float64 `json:"zi"`
	} `json:"params"`
	Output struct {
		Y  []float64 `json:"y"`
		Zf []float64 `json:"zf"`
	} `json:"output"`
}

// filterStateSciPySuite checks dsp.FilterState against scipy.signal.lfilter with zi,
// including the final state
// measured: worst 2.3e-16 (1 ULP) across all six filters, including a[0] != 1,
// which confirms SciPy treats zi as the state of the normalised filter, as apply
// does; set to match Filter
var filterStateSciPySuite = newSuiteOf("FilterState", reference.SourceSciPy, "lfilter_state",
	reference.Tolerance{Rel: 1e-14, Scaled: 1e-14},
	filterStateOutputs, isWellConditioned)

// filterStateMATLABSuite checks dsp.FilterState against MATLAB's [y, zf] = filter(b, a, x, zi)
// measured: worst 3.0e-16 (1 ULP), confirming MATLAB's zi is the same transposed
// direct form II state as SciPy's; set to match Filter
var filterStateMATLABSuite = newSuiteOf("FilterState", reference.SourceMATLAB, "lfilter_state",
	reference.Tolerance{Rel: 1e-14, Scaled: 1e-14},
	filterStateOutputs, isWellConditioned)

// filterStateOutputs pairs dsp.FilterState's output and final state with a vector's
func filterStateOutputs(v lfilterStateVector) (outputs []Output, err error) {
	y, zf, err := dsp.FilterState(v.Params.B, v.Params.A, v.Params.X, v.Params.Zi)
	if err != nil {
		return nil, err
	}
	scale := maxAbs(v.Params.X)
	outputs = []Output{
		{Name: "y", Got: y, Want: v.Output.Y, Scale: scale},
		{Name: "zf", Got: zf, Want: v.Output.Zf, Scale: scale},
	}
	return
}

type lfilterZiVector struct {
	reference.Meta
	Params struct {
		B []float64 `json:"b"`
		A []float64 `json:"a"`
	} `json:"params"`
	Output struct {
		Zi []float64 `json:"zi"`
	} `json:"output"`
}

// filterZiSciPySuite checks dsp.FilterZi against scipy.signal.lfilter_zi
// MATLAB has no equivalent function
// measured: exact on amd64, and worst 3.3e-16 (2 ULP) on arm64, on butter2_lp0.001;
// set to match Filter
var filterZiSciPySuite = newSuite("FilterZi", reference.SourceSciPy, "lfilter_zi",
	reference.Tolerance{Rel: 1e-14, Scaled: 1e-14},
	func(v lfilterZiVector) (outputs []Output, err error) {
		zi, err := dsp.FilterZi(v.Params.B, v.Params.A)
		if err != nil {
			return nil, err
		}
		outputs = []Output{{Name: "zi", Got: zi, Want: v.Output.Zi}}
		return
	})
