package conformance

import (
	"fmt"
	"math"
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

type savGolVector struct {
	reference.Meta
	Params struct {
		X         reference.Floats `json:"x"` // may hold NaN
		N         int              `json:"n"`
		PolyOrder int              `json:"polyorder"`
		Deriv     int              `json:"deriv"`
		Delta     float64          `json:"delta"`
	} `json:"params"`
	Output struct {
		Y reference.Floats `json:"y"`
	} `json:"output"`
}

// SavGolFilter, measured on both architectures. The SciPy vectors are the
// exact filter: rational least-squares weights, each output rounded once,
// rather than scipy.signal.savgol_filter's output, which fits raw powers of
// the positions by lstsq and fails at n = 101, degree 10, wrong by 21%.
// SavGolFilter is worst at 3.6e-16 (2 ULP of the input's scale over
// delta^deriv), on chirp_n11_p4_d2, and at 2.3e-16 on long101_n101_p10.
//
// MATLAB's sgolayfilt only smooths, so its vectors are the smoothing ones,
// worst at 8.7e-16 on long51_n51_p6. At n = 101, degree 10, it loses digits,
// 2.5e-14 off, so long101_n101_p10 is in a suite of its own.
var (
	savGolSciPySuite  = newSuite("SavGolFilter", reference.SourceSciPy, "savgol_filter", tolerance(1e-14), savGolOutputs)
	savGolMATLABSuite = newSuiteOf("SavGolFilter", reference.SourceMATLAB, "savgol_filter", tolerance(1e-14), savGolOutputs,
		func(file string) bool { return !longSavGol(file) })
	savGolLongMATLABSuite = newSuiteOf("SavGolFilter (n = 101)", reference.SourceMATLAB, "savgol_filter", tolerance(1e-12), savGolOutputs,
		longSavGol)
)

// longSavGol accepts the vectors of a 101-sample window and degree 10, where
// MATLAB's sgolayfilt loses digits
func longSavGol(file string) bool {
	return strings.HasPrefix(file, "long101_")
}

// savGolOutputs pairs dsp.SavGolFilter's output with a vector's, measured
// against the scale of the input over delta^deriv, which an output's rounding
// grows with even where it's ~0, as a derivative of a slow signal is
func savGolOutputs(v savGolVector) (outputs []Output, err error) {
	y, err := dsp.SavGolFilter(v.Params.X, v.Params.N, v.Params.PolyOrder, dsp.SavGolDeriv(v.Params.Deriv, v.Params.Delta))
	if err != nil {
		return nil, err
	}
	scale := finiteMaxAbs(v.Params.X) / math.Pow(v.Params.Delta, float64(v.Params.Deriv))
	outputs = []Output{{Name: "y", Got: y, Want: v.Output.Y, Scale: scale}}
	return
}

// finiteMaxAbs is the largest magnitude in v, leaving out NaN
func finiteMaxAbs(v []float64) (m float64) {
	for _, x := range v {
		if !math.IsNaN(x) {
			m = math.Max(m, math.Abs(x))
		}
	}
	return
}
