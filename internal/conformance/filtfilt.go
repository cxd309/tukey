package conformance

import (
	"fmt"

	"github.com/cxd309/tukey/dsp"
	"github.com/cxd309/tukey/internal/reference"
)

type filtfiltVector struct {
	reference.Meta
	Params struct {
		B       []float64 `json:"b"`
		A       []float64 `json:"a"`
		X       []float64 `json:"x"`
		PadType *string   `json:"padtype"` // null: no padding
		PadLen  *int      `json:"padlen"`  // null: SciPy's default
	} `json:"params"`
	Output struct {
		Y []float64 `json:"y"`
	} `json:"output"`
}

// filtfiltSciPySuite checks dsp.FiltFilt against scipy.signal.filtfilt
// measured: worst 7.7e-16 (3 ULP); set ~10x above
var filtfiltSciPySuite = newSuiteOf("FiltFilt", reference.SourceSciPy, "filtfilt",
	reference.Tolerance{Rel: 1e-14, Scaled: 1e-14},
	filtfiltOutputs(), isWellConditioned)

// filtfiltMATLABSuite checks dsp.FiltFilt with MATLABPadLen against MATLAB's filtfilt
// MATLAB has one edge handling, so its vectors are SciPy's default-padding cases only
// measured: worst 7.7e-16 (3 ULP) on all 37 vectors, which verifies MATLABPadLen;
// set to match the SciPy suite
var filtfiltMATLABSuite = newSuiteOf("FiltFilt", reference.SourceMATLAB, "filtfilt",
	reference.Tolerance{Rel: 1e-14, Scaled: 1e-14},
	filtfiltOutputs(dsp.MATLABPadLen()), isWellConditioned)

// filtfiltOutputs pairs dsp.FiltFilt's output with a vector's,
// with the vector's padding plus any extra options
func filtfiltOutputs(extra ...dsp.FiltOption) (check func(v filtfiltVector) ([]Output, error)) {
	check = func(v filtfiltVector) (outputs []Output, err error) {
		opts, err := filtfiltOptions(v.Params.PadType, v.Params.PadLen)
		if err != nil {
			return nil, err
		}
		y, err := dsp.FiltFilt(v.Params.B, v.Params.A, v.Params.X, append(opts, extra...)...)
		if err != nil {
			return nil, err
		}
		outputs = []Output{{Name: "y", Got: y, Want: v.Output.Y, Scale: maxAbs(v.Params.X)}}
		return
	}
	return
}

// filtfiltOptions translates a vector's SciPy padtype and padlen into FiltOptions
func filtfiltOptions(padtype *string, padlen *int) (opts []dsp.FiltOption, err error) {
	padding := dsp.PaddingNone
	if padtype != nil {
		switch *padtype {
		case "odd":
			padding = dsp.PaddingOdd
		case "even":
			padding = dsp.PaddingEven
		case "constant":
			padding = dsp.PaddingConstant
		default:
			return nil, fmt.Errorf("unknown padtype %q", *padtype)
		}
	}
	opts = append(opts, dsp.PadType(padding))
	if padlen != nil {
		opts = append(opts, dsp.PadLen(*padlen))
	}
	return
}
