package conformance

import (
	"fmt"

	"github.com/cxd309/tukey/dsp"
	"github.com/cxd309/tukey/internal/reference"
)

type detrendVector struct {
	reference.Meta
	Params struct {
		X    []float64 `json:"x"`
		Type string    `json:"type"`
	} `json:"params"`
	Output struct {
		Y []float64 `json:"y"`
	} `json:"output"`
}

// detrend measured: worst 5.6e-16 (2 ULP) against SciPy and 1.2e-15 (5 ULP)
// against MATLAB, on detrending a constant, where the fitted line cancels the
// signal almost exactly. Set ~10x above the worst.
var (
	detrendSciPySuite  = newSuite("Detrend", reference.SourceSciPy, "detrend", tolerance(1e-14), detrendOutputs)
	detrendMATLABSuite = newSuite("Detrend", reference.SourceMATLAB, "detrend", tolerance(1e-13), detrendOutputs)
)

// detrendOutputs pairs dsp.Detrend's output with a vector's. It's scaled to the
// input, as detrending a nearly straight signal cancels to ~0.
func detrendOutputs(v detrendVector) (outputs []Output, err error) {
	mode := dsp.DetrendLinear
	switch v.Params.Type {
	case "linear":
	case "constant":
		mode = dsp.DetrendConstant
	default:
		return nil, fmt.Errorf("unknown detrend type %q", v.Params.Type)
	}
	y, err := dsp.Detrend(v.Params.X, mode)
	if err != nil {
		return nil, err
	}
	outputs = []Output{{Name: "y", Got: y, Want: v.Output.Y, Scale: maxAbs(v.Params.X)}}
	return
}

type interp1Vector struct {
	reference.Meta
	Params struct {
		X  []float64 `json:"x"`
		Y  []float64 `json:"y"`
		XI []float64 `json:"xi"`
	} `json:"params"`
	Output struct {
		YI reference.Floats `json:"yi"` // NaN outside x's range
	} `json:"output"`
}

// interp1 measured: worst 1.3e-16 against SciPy and 2.2e-16 against MATLAB
// (1 ULP), with NaN outside x's range in the same places as both. Set ~10x above.
var (
	interp1SciPySuite  = newSuite("Interp1", reference.SourceSciPy, "interp1", tolerance(1e-14), interp1Outputs)
	interp1MATLABSuite = newSuite("Interp1", reference.SourceMATLAB, "interp1", tolerance(1e-14), interp1Outputs)
)

// interp1Outputs pairs dsp.Interp1's output with a vector's.
func interp1Outputs(v interp1Vector) (outputs []Output, err error) {
	yi, err := dsp.Interp1(v.Params.X, v.Params.Y, v.Params.XI)
	if err != nil {
		return nil, err
	}
	outputs = []Output{{Name: "yi", Got: yi, Want: v.Output.YI}}
	return
}

type fillMissingVector struct {
	reference.Meta
	Params struct {
		X      reference.Floats `json:"x"`
		Method string           `json:"method"`
		Value  *float64         `json:"value"` // the constant, or null for linear
	} `json:"params"`
	Output struct {
		Y reference.Floats `json:"y"`
	} `json:"output"`
}

// fillmissing measured: exact against both, NaN included, on both architectures.
// The SciPy side is MATLAB's documented definition, which MATLAB's own vectors
// confirm. Exact here is a measurement, not a guarantee (the same line formula,
// rounded the same way), so the tolerance stays at 1e-14.
var (
	fillMissingSciPySuite  = newSuite("FillMissing", reference.SourceSciPy, "fillmissing", tolerance(1e-14), fillMissingOutputs)
	fillMissingMATLABSuite = newSuite("FillMissing", reference.SourceMATLAB, "fillmissing", tolerance(1e-14), fillMissingOutputs)
)

// fillMissingOutputs pairs dsp.FillMissing's output with a vector's. It sets no
// Scale, because the inputs hold NaN and maxAbs would make it NaN.
func fillMissingOutputs(v fillMissingVector) (outputs []Output, err error) {
	var fill dsp.Fill
	switch v.Params.Method {
	case "linear":
		fill = dsp.FillLinear()
	case "constant":
		if v.Params.Value == nil {
			return nil, fmt.Errorf("constant fill without a value")
		}
		fill = dsp.FillConstant(*v.Params.Value)
	default:
		return nil, fmt.Errorf("unknown fill method %q", v.Params.Method)
	}
	y := dsp.FillMissing(v.Params.X, fill)
	outputs = []Output{{Name: "y", Got: y, Want: v.Output.Y}}
	return
}

type rmsVector struct {
	reference.Meta
	Params struct {
		X []float64 `json:"x"`
	} `json:"params"`
	Output struct {
		RMS []float64 `json:"rms"`
	} `json:"output"`
}

// rms measured: worst 1.4e-16 (1 ULP) against SciPy and 1.2e-16 (1 ULP) against MATLAB,
// on both architectures, with the squares summed with compensation. Set ~10x above.
var (
	rmsSciPySuite  = newSuite("RMS", reference.SourceSciPy, "rms", tolerance(1e-14), rmsOutputs)
	rmsMATLABSuite = newSuite("RMS", reference.SourceMATLAB, "rms", tolerance(1e-14), rmsOutputs)
)

// rmsOutputs pairs dsp.RMS with a vector's.
func rmsOutputs(v rmsVector) (outputs []Output, err error) {
	rms, err := dsp.RMS(v.Params.X)
	if err != nil {
		return nil, err
	}
	outputs = []Output{{Name: "rms", Got: []float64{rms}, Want: v.Output.RMS}}
	return
}
