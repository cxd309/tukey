package conformance

import (
	"github.com/cxd309/tukey/dsp"
	"github.com/cxd309/tukey/internal/reference"
)

type hilbertVector struct {
	reference.Meta
	Params struct {
		X []float64 `json:"x"`
	} `json:"params"`
	Output struct {
		Analytic reference.Complex `json:"analytic"`
	} `json:"output"`
}

// hilbert measured: worst 5.7e-16 (3 ULP) against SciPy and 5.8e-16 (3 ULP)
// against MATLAB, on the longest signal, on both architectures. Set ~10x above.
var (
	hilbertSciPySuite  = newSuite("Hilbert", reference.SourceSciPy, "hilbert", tolerance(1e-14), hilbertOutputs)
	hilbertMATLABSuite = newSuite("Hilbert", reference.SourceMATLAB, "hilbert", tolerance(1e-14), hilbertOutputs)
)

// hilbertOutputs pairs dsp.Hilbert's analytic signal with a vector's.
func hilbertOutputs(v hilbertVector) (outputs []Output, err error) {
	analytic, err := dsp.Hilbert(v.Params.X)
	if err != nil {
		return nil, err
	}
	want, err := v.Output.Analytic.Values()
	if err != nil {
		return nil, err
	}
	outputs = []Output{{Name: "analytic", GotComplex: analytic, WantComplex: want}}
	return
}
