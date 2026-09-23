package conformance

import (
	"github.com/cxd309/godsp/dsp"
	"github.com/cxd309/godsp/internal/reference"
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

// filterSuite checks dsp.Filter against scipy.signal.lfilter
// measured: FIR and simple-coefficient cases are bit-identical; Butterworth cases
// agree within ~3 ULPs (worst 6.6e-16), differing only in the order of additions;
// set ~10x above the worst for headroom across platforms
var filterSuite = newSuite("Filter", "filter",
	reference.Tolerance{Rel: 1e-14, Abs: 1e-14},
	func(v filterVector) (outputs []Output, err error) {
		y, err := dsp.Filter(v.Params.B, v.Params.A, v.Params.X)
		if err != nil {
			return nil, err
		}
		outputs = []Output{{"y", y, v.Output.Y}}
		return
	})
