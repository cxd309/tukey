package conformance

import (
	"github.com/cxd309/tukey/dsp"
	"github.com/cxd309/tukey/internal/reference"
)

type sosfiltVector struct {
	reference.Meta
	Params struct {
		SOS [][6]float64 `json:"sos"`
		X   []float64    `json:"x"`
		Zi  [][2]float64 `json:"zi"` // null: start from rest
	} `json:"params"`
	Output struct {
		Y  []float64    `json:"y"`
		Zf [][2]float64 `json:"zf"` // null when zi was
	} `json:"output"`
}

// sosFilterSciPySuite checks dsp.SOSFilter against scipy.signal.sosfilt, and
// dsp.SOSFilterState against sosfilt with zi, including the final state
// measured: worst 4.1e-14 (~183 ULP) on arm64 and 1.3e-15 on amd64, both on
// butter8_lp0.01, whose poles sit very close to z = 1, so its feedback amplifies
// single-rounding differences. The arm64 gap is fused multiply-add: Go fuses on
// arm64 and NumPy doesn't. Other filters, including the narrow order-20 bandpass,
// agree within a few ULP. Set ~10x above the worst.
var sosFilterSciPySuite = newSuite("SOSFilter", reference.SourceSciPy, "sosfilt",
	reference.Tolerance{Rel: 1e-12, Scaled: 1e-12},
	sosFilterOutputs)

// sosFilterMATLABSuite checks dsp.SOSFilter against MATLAB's sosfilt
// sosfilt takes no initial state, so the with-state cases have no MATLAB vector
// measured: worst 9.0e-15 (41 ULP) on arm64 and 4.5e-14 on amd64, both on
// butter8_lp0.01. That's the mirror image of the SciPy suite: tukey agrees best
// with MATLAB when it fuses multiply-adds, so MATLAB evidently does too.
// Set ~10x above the worst, the same as the SciPy suite.
var sosFilterMATLABSuite = newSuite("SOSFilter", reference.SourceMATLAB, "sosfilt",
	reference.Tolerance{Rel: 1e-12, Scaled: 1e-12},
	sosFilterOutputs)

// sosFilterOutputs pairs dsp.SOSFilter's output with a vector's,
// or dsp.SOSFilterState's output and final state when the vector has a starting state
func sosFilterOutputs(v sosfiltVector) (outputs []Output, err error) {
	scale := maxAbs(v.Params.X)
	if v.Params.Zi == nil {
		y, err := dsp.SOSFilter(v.Params.SOS, v.Params.X)
		if err != nil {
			return nil, err
		}
		outputs = []Output{{Name: "y", Got: y, Want: v.Output.Y, Scale: scale}}
		return outputs, nil
	}
	y, zf, err := dsp.SOSFilterState(v.Params.SOS, v.Params.X, v.Params.Zi)
	if err != nil {
		return nil, err
	}
	outputs = []Output{
		{Name: "y", Got: y, Want: v.Output.Y, Scale: scale},
		{Name: "zf", Got: flattenPairs(zf), Want: flattenPairs(v.Output.Zf), Scale: scale},
	}
	return
}

type sosfiltZiVector struct {
	reference.Meta
	Params struct {
		SOS [][6]float64 `json:"sos"`
	} `json:"params"`
	Output struct {
		Zi [][2]float64 `json:"zi"`
	} `json:"output"`
}

// sosFilterZiSciPySuite checks dsp.SOSFilterZi against scipy.signal.sosfilt_zi
// MATLAB has no equivalent function
// measured: exact or under 1 ULP on all 4 filters (worst 9.8e-17); set to match Filter
var sosFilterZiSciPySuite = newSuite("SOSFilterZi", reference.SourceSciPy, "sosfilt_zi",
	reference.Tolerance{Rel: 1e-14, Scaled: 1e-14},
	func(v sosfiltZiVector) (outputs []Output, err error) {
		zi, err := dsp.SOSFilterZi(v.Params.SOS)
		if err != nil {
			return nil, err
		}
		outputs = []Output{{Name: "zi", Got: flattenPairs(zi), Want: flattenPairs(v.Output.Zi)}}
		return
	})

type sosfiltfiltVector struct {
	reference.Meta
	Params struct {
		SOS     [][6]float64 `json:"sos"`
		X       []float64    `json:"x"`
		PadType *string      `json:"padtype"` // null: no padding
		PadLen  *int         `json:"padlen"`  // null: SciPy's default
	} `json:"params"`
	Output struct {
		Y []float64 `json:"y"`
	} `json:"output"`
}

// sosFiltFiltSciPySuite checks dsp.SOSFiltFilt against scipy.signal.sosfiltfilt
// measured: worst 1.7e-13 (~763 ULP) on butter8_lp0.01: SOSFilter's sensitivity on
// that filter, over two passes of a padded signal; other filters within a few ULP;
// set ~10x above the worst
//
// not yet checked against MATLAB: its filtfilt(sos, g, x) runs forward-backward per
// section, each padded by 3 * the total filter order, which SOSFiltFilt doesn't
// offer; testdata/matlab/sosfiltfilt is waiting for the MATLAB compatibility
// option on the roadmap
var sosFiltFiltSciPySuite = newSuite("SOSFiltFilt", reference.SourceSciPy, "sosfiltfilt",
	reference.Tolerance{Rel: 1e-11, Scaled: 1e-11},
	func(v sosfiltfiltVector) (outputs []Output, err error) {
		opts, err := filtfiltOptions(v.Params.PadType, v.Params.PadLen)
		if err != nil {
			return nil, err
		}
		y, err := dsp.SOSFiltFilt(v.Params.SOS, v.Params.X, opts...)
		if err != nil {
			return nil, err
		}
		outputs = []Output{{Name: "y", Got: y, Want: v.Output.Y, Scale: maxAbs(v.Params.X)}}
		return
	})

// flattenPairs lays per-section state pairs out as one slice, for comparison
func flattenPairs(pairs [][2]float64) (flat []float64) {
	for _, p := range pairs {
		flat = append(flat, p[0], p[1])
	}
	return
}
