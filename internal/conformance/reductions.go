package conformance

import (
	"fmt"
	"math"
	"strings"

	"github.com/cxd309/tukey/dsp"
	"github.com/cxd309/tukey/internal/reference"
)

type sumVector struct {
	reference.Meta
	Params struct {
		X []float64 `json:"x"`
	} `json:"params"`
	Output struct {
		Sum  []float64 `json:"sum"`
		Mean []float64 `json:"mean"`
		S    []float64 `json:"s"` // cumsum
	} `json:"output"`
}

// sum, mean and cumsum measured on both architectures: exact against SciPy, whose side is
// summed exactly, so it's measured against each value itself and the tolerance
// stays at 1e-14. MATLAB sums plainly, 5.1e-7 from exact on mixed_magnitudes
// relative to the total, so it's measured against the size of what's summed,
// sum(|x|), the scale it rounds at. Against that, sum and mean are worst at
// 9.0e-14 (407 ULP) and cumsum at 3.6e-13 (1629 ULP), all on tenths_long, as a
// plain running sum's error grows with its length. Set ~10x above the worst.
var (
	sumSciPySuite     = newSuite("Sum", reference.SourceSciPy, "sum", tolerance(1e-14), sumOutputs(false))
	sumMATLABSuite    = newSuite("Sum", reference.SourceMATLAB, "sum", tolerance(1e-12), sumOutputs(true))
	meanSciPySuite    = newSuite("Mean", reference.SourceSciPy, "mean", tolerance(1e-14), meanOutputs(false))
	meanMATLABSuite   = newSuite("Mean", reference.SourceMATLAB, "mean", tolerance(1e-12), meanOutputs(true))
	cumSumSciPySuite  = newSuite("CumSum", reference.SourceSciPy, "cumsum", tolerance(1e-14), cumSumOutputs(false))
	cumSumMATLABSuite = newSuite("CumSum", reference.SourceMATLAB, "cumsum", tolerance(1e-11), cumSumOutputs(true))
)

// sumOutputs pairs dsp.Sum with a vector's. scaled measures against sum(|x|).
func sumOutputs(scaled bool) func(v sumVector) ([]Output, error) {
	return func(v sumVector) (outputs []Output, err error) {
		outputs = []Output{{Name: "sum", Got: []float64{dsp.Sum(v.Params.X)}, Want: v.Output.Sum, Scale: sumScale(v.Params.X, scaled)}}
		return
	}
}

// meanOutputs pairs dsp.Mean with a vector's. scaled measures against mean(|x|).
func meanOutputs(scaled bool) func(v sumVector) ([]Output, error) {
	return func(v sumVector) (outputs []Output, err error) {
		m, err := dsp.Mean(v.Params.X)
		if err != nil {
			return nil, err
		}
		scale := sumScale(v.Params.X, scaled) / float64(len(v.Params.X))
		outputs = []Output{{Name: "mean", Got: []float64{m}, Want: v.Output.Mean, Scale: scale}}
		return
	}
}

// cumSumOutputs pairs dsp.CumSum with a vector's. scaled measures against sum(|x|).
func cumSumOutputs(scaled bool) func(v sumVector) ([]Output, error) {
	return func(v sumVector) (outputs []Output, err error) {
		outputs = []Output{{Name: "s", Got: dsp.CumSum(v.Params.X), Want: v.Output.S, Scale: sumScale(v.Params.X, scaled)}}
		return
	}
}

type diffVector struct {
	reference.Meta
	Params struct {
		X []float64 `json:"x"`
	} `json:"params"`
	Output struct {
		D []float64 `json:"d"`
	} `json:"output"`
}

// diff measured on both architectures: exact against both. Each difference is one
// subtraction, the same in all three, so exact is expected, but the tolerance
// stays at 1e-14.
var (
	diffSciPySuite  = newSuite("Diff", reference.SourceSciPy, "diff", tolerance(1e-14), diffOutputs)
	diffMATLABSuite = newSuite("Diff", reference.SourceMATLAB, "diff", tolerance(1e-14), diffOutputs)
)

// diffOutputs pairs dsp.Diff with a vector's.
func diffOutputs(v diffVector) (outputs []Output, err error) {
	outputs = []Output{{Name: "d", Got: dsp.Diff(v.Params.X), Want: v.Output.D}}
	return
}

type gradVector struct {
	reference.Meta
	Params struct {
		X  []float64 `json:"x"`
		DX float64   `json:"dx"`
	} `json:"params"`
	Output struct {
		G []float64 `json:"g"`
	} `json:"output"`
}

// gradient measured on both architectures: exact against SciPy, as Grad's arithmetic is
// NumPy's, so the tolerance stays at 1e-14. Worst 2.1e-14 (96 ULP) against
// MATLAB on sine, whose dx is 0.01, as MATLAB rounds its differences
// differently. Set ~10x above the worst.
var (
	gradSciPySuite  = newSuite("Grad", reference.SourceSciPy, "gradient", tolerance(1e-14), gradOutputs)
	gradMATLABSuite = newSuite("Grad", reference.SourceMATLAB, "gradient", tolerance(1e-12), gradOutputs)
)

// gradOutputs pairs dsp.Grad with a vector's.
func gradOutputs(v gradVector) (outputs []Output, err error) {
	g, err := dsp.Grad(v.Params.X, v.Params.DX)
	if err != nil {
		return nil, err
	}
	outputs = []Output{{Name: "g", Got: g, Want: v.Output.G}}
	return
}

type trapzVector struct {
	reference.Meta
	Params struct {
		Y  []float64 `json:"y"`
		DX float64   `json:"dx"`
	} `json:"params"`
	Output struct {
		Integral []float64 `json:"integral"`
	} `json:"output"`
}

// trapezoid and cumulative_trapezoid measured on both architectures: exact against SciPy,
// whose side is summed exactly as for sum, so the tolerance stays at 1e-14.
// MATLAB is measured against dx·sum(|y|): worst 9.0e-14 (407 ULP) for trapz and
// 3.6e-13 (1629 ULP) for cumtrapz, both on tenths_long, as for sum and cumsum.
// Set ~10x above the worst.
var (
	trapzSciPySuite     = newSuite("Trapz", reference.SourceSciPy, "trapezoid", tolerance(1e-14), trapzOutputs(false))
	trapzMATLABSuite    = newSuite("Trapz", reference.SourceMATLAB, "trapezoid", tolerance(1e-12), trapzOutputs(true))
	cumTrapzSciPySuite  = newSuite("CumTrapz", reference.SourceSciPy, "cumulative_trapezoid", tolerance(1e-14), cumTrapzOutputs(false))
	cumTrapzMATLABSuite = newSuite("CumTrapz", reference.SourceMATLAB, "cumulative_trapezoid", tolerance(1e-11), cumTrapzOutputs(true))
)

// trapzOutputs pairs dsp.Trapz with a vector's. scaled measures against dx·sum(|y|).
func trapzOutputs(scaled bool) func(v trapzVector) ([]Output, error) {
	return func(v trapzVector) (outputs []Output, err error) {
		integral, err := dsp.Trapz(v.Params.Y, v.Params.DX)
		if err != nil {
			return nil, err
		}
		scale := v.Params.DX * sumScale(v.Params.Y, scaled)
		outputs = []Output{{Name: "integral", Got: []float64{integral}, Want: v.Output.Integral, Scale: scale}}
		return
	}
}

// cumTrapzOutputs pairs dsp.CumTrapz with a vector's. scaled measures against dx·sum(|y|).
func cumTrapzOutputs(scaled bool) func(v trapzVector) ([]Output, error) {
	return func(v trapzVector) (outputs []Output, err error) {
		integral, err := dsp.CumTrapz(v.Params.Y, v.Params.DX)
		if err != nil {
			return nil, err
		}
		scale := v.Params.DX * sumScale(v.Params.Y, scaled)
		outputs = []Output{{Name: "integral", Got: integral, Want: v.Output.Integral, Scale: scale}}
		return
	}
}

type percentileVector struct {
	reference.Meta
	Params struct {
		X      reference.Floats `json:"x"` // may hold NaN
		P      []float64        `json:"p"`
		Method string           `json:"method"`
	} `json:"params"`
	Output struct {
		Values reference.Floats `json:"values"`
	} `json:"output"`
}

// percentile measured on both architectures: exact against SciPy, as the interpolation is
// NumPy's, and within 1 ULP of MATLAB (worst 4.0e-17 on noise101_hazen), so the
// tolerance stays at 1e-14. MATLAB's prctile is NumPy's "hazen", and ignores NaN
// where tukey and NumPy give NaN, so the NaN case isn't compared against MATLAB.
var (
	percentileSciPySuite  = newSuite("Percentile", reference.SourceSciPy, "percentile", tolerance(1e-14), percentileOutputs)
	percentileMATLABSuite = newSuiteOf("Percentile", reference.SourceMATLAB, "percentile", tolerance(1e-14), percentileOutputs,
		func(file string) bool { return !strings.HasPrefix(file, "nan_") })
)

// percentileOutputs pairs dsp.Percentile at each of a vector's p with its values.
func percentileOutputs(v percentileVector) (outputs []Output, err error) {
	var opts []dsp.PercentileOption
	switch v.Params.Method {
	case "linear":
	case "hazen":
		opts = append(opts, dsp.PercentileBy(dsp.PercentileHazen))
	default:
		return nil, fmt.Errorf("unknown percentile method %q", v.Params.Method)
	}
	got := make([]float64, len(v.Params.P))
	for i, p := range v.Params.P {
		if got[i], err = dsp.Percentile(v.Params.X, p, opts...); err != nil {
			return nil, err
		}
	}
	outputs = []Output{{Name: "values", Got: got, Want: v.Output.Values}}
	return
}

// sumScale is sum(|x|), the size of what's summed and so the scale a plain sum
// rounds at, if scaled, or 0 to measure against the values themselves.
func sumScale(x []float64, scaled bool) (scale float64) {
	if !scaled {
		return 0
	}
	for _, v := range x {
		scale += math.Abs(v)
	}
	return
}
