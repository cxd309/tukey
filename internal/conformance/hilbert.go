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

// The hilbert suites.
// Measured: worst 3 ULP, 5.7e-16 against SciPy and 5.8e-16 against MATLAB, on the longest signal, on both architectures.
// Tolerance: about 10 times the worst.
var (
	hilbertSciPySuite  = newSuite("Hilbert", reference.SourceSciPy, "hilbert", tolerance(1e-14), hilbertOutputs)
	hilbertMATLABSuite = newSuite("Hilbert", reference.SourceMATLAB, "hilbert", tolerance(1e-14), hilbertOutputs)
)

// hilbertOutputs pairs dsp.Hilbert's analytic signal with a vector's.
func hilbertOutputs(v hilbertVector) (outputs []Output, err error) {
	analytic := dsp.Hilbert(v.Params.X)
	want, err := v.Output.Analytic.Values()
	if err != nil {
		return nil, err
	}
	outputs = []Output{{Name: "analytic", GotComplex: analytic, WantComplex: want}}
	return
}
