package conformance

import (
	"github.com/cxd309/tukey/dsp"
	"github.com/cxd309/tukey/internal/reference"
)

type sosFiltFiltBySectionVector struct {
	reference.Meta
	Params struct {
		SOS          [][6]float64 `json:"sos"`
		X            []float64    `json:"x"`
		PadType      *string      `json:"padtype"`       // or null for no padding
		PadLen       *int         `json:"padlen"`        // or null for SciPy's default, or MATLAB's
		MATLABPadLen bool         `json:"matlab_padlen"` // MATLAB's padding length
	} `json:"params"`
	Output struct {
		Y []float64 `json:"y"`
	} `json:"output"`
}

type filtfiltMATLABVector struct {
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

// SOSFiltFilt with BySection, and FiltFilt with MATLABPadLen, measured on both architectures.
// MATLAB's filtfilt filters each section forward and backward in turn, padded by 3 × the order,
// which it counts to each polynomial's last nonzero coefficient.
// SciPy has no equivalent, so its side is that definition built from scipy.signal.filtfilt, one section at a time,
// which MATLAB's own vectors agree with to rounding.
//
// By section, against SciPy's side, it's worst at 7.7e-13 on butter8_lp0.01, the difficult filter, same as SOSFilter,
// over two passes of each of its four sections.
// The rest, the order-count filters included, are within 2.1e-14.
// Against MATLAB it's worst at 7.0e-13, on the same filter.
// MATLAB's vectors are the cases with MATLAB's padding, all 25 of testdata/matlab/sosfiltfilt among them, identical,
// so that category isn't compared again.
//
// FiltFilt with MATLABPadLen, on (b, a) with trailing zeros and a pure gain, whose edge is MATLAB's minimum of 1,
// is worst at 2 ULP, 4.3e-16, against both.
var (
	sosFiltFiltBySectionSciPySuite = newSuite("SOSFiltFilt + BySection", reference.SourceSciPy, "sosfiltfilt_by_section",
		tolerance(1e-11), sosFiltFiltBySectionOutputs)
	sosFiltFiltBySectionMATLABSuite = newSuite("SOSFiltFilt + BySection", reference.SourceMATLAB, "sosfiltfilt_by_section",
		tolerance(1e-11), sosFiltFiltBySectionOutputs)
	filtfiltMATLABPadLenSciPySuite = newSuite("FiltFilt + MATLABPadLen", reference.SourceSciPy, "filtfilt_matlab",
		tolerance(1e-14), filtfiltMATLABPadLenOutputs)
	filtfiltMATLABPadLenMATLABSuite = newSuite("FiltFilt + MATLABPadLen", reference.SourceMATLAB, "filtfilt_matlab",
		tolerance(1e-14), filtfiltMATLABPadLenOutputs)
)

// sosFiltFiltBySectionOutputs pairs dsp.SOSFiltFilt, with BySection and a vector's padding, with its output.
func sosFiltFiltBySectionOutputs(v sosFiltFiltBySectionVector) (outputs []Output, err error) {
	opts, err := filtfiltOptions(v.Params.PadType, v.Params.PadLen)
	if err != nil {
		return nil, err
	}
	opts = append(opts, dsp.BySection())
	if v.Params.MATLABPadLen {
		opts = append(opts, dsp.MATLABPadLen())
	}
	y, err := dsp.SOSFiltFilt(v.Params.SOS, v.Params.X, opts...)
	if err != nil {
		return nil, err
	}
	outputs = []Output{{Name: "y", Got: y, Want: v.Output.Y, Scale: maxAbs(v.Params.X)}}
	return
}

// filtfiltMATLABPadLenOutputs pairs dsp.FiltFilt, with MATLABPadLen, with a vector's output.
func filtfiltMATLABPadLenOutputs(v filtfiltMATLABVector) (outputs []Output, err error) {
	y, err := dsp.FiltFilt(v.Params.B, v.Params.A, v.Params.X, dsp.MATLABPadLen())
	if err != nil {
		return nil, err
	}
	outputs = []Output{{Name: "y", Got: y, Want: v.Output.Y, Scale: maxAbs(v.Params.X)}}
	return
}
