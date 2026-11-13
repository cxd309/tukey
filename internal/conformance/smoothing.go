package conformance

import (
	"fmt"
	"strings"

	"github.com/cxd309/tukey/dsp"
	"github.com/cxd309/tukey/internal/reference"
)

type medFiltVector struct {
	reference.Meta
	Params struct {
		X    reference.Floats `json:"x"` // may hold NaN
		N    int              `json:"n"`
		Edge string           `json:"edge"`
	} `json:"params"`
	Output struct {
		Y reference.Floats `json:"y"`
	} `json:"output"`
}

// MedFilt, measured on both architectures. The SciPy vectors are the median
// filter from its definition, with an even-count window's midpoint exact,
// and match scipy.signal.medfilt where it applies: zero-padded, without NaN.
// MedFilt is exact against all of them, as a median only selects, and its
// midpoint is rounded once.
//
// MATLAB's medfilt1 agrees except at some even-count midpoints near the
// truncated ends, 1 ULP off, worst at 2.8e-17 on chirp_n11_truncate. Its
// midpoint of two subnormals rounds twice, as a/2 + b/2 would:
// subnormal_n3_truncate and subnormal_n5_truncate are 1 subnormal ULP, half
// their own value, off, and aren't compared.
var (
	medFiltSciPySuite  = newSuite("MedFilt", reference.SourceSciPy, "medfilt", tolerance(1e-14), medFiltOutputs)
	medFiltMATLABSuite = newSuiteOf("MedFilt", reference.SourceMATLAB, "medfilt", tolerance(1e-14), medFiltOutputs,
		func(file string) bool { return !strings.HasPrefix(file, "subnormal_") })
)

// medFiltOutputs pairs dsp.MedFilt's output with a vector's. It sets no
// Scale, because the inputs hold NaN and maxAbs would make it NaN.
func medFiltOutputs(v medFiltVector) (outputs []Output, err error) {
	var edge dsp.MedFiltEdge
	switch v.Params.Edge {
	case "truncate":
		edge = dsp.MedFiltTruncate
	case "zeropad":
		edge = dsp.MedFiltZeroPad
	default:
		return nil, fmt.Errorf("unknown medfilt edge %q", v.Params.Edge)
	}
	y, err := dsp.MedFilt(v.Params.X, v.Params.N, dsp.MedFiltEdges(edge))
	if err != nil {
		return nil, err
	}
	outputs = []Output{{Name: "y", Got: y, Want: v.Output.Y}}
	return
}
