package conformance

import (
	"encoding/json"
	"fmt"

	"github.com/cxd309/tukey/dsp"
	"github.com/cxd309/tukey/internal/reference"
)

type besselVector struct {
	reference.Meta
	Params struct {
		Order int             `json:"order"`
		Wn    json.RawMessage `json:"wn"`
		BType string          `json:"btype"`
		Norm  string          `json:"norm"` // SciPy's: "phase", "delay" or "mag"
	} `json:"params"`
	Output struct {
		B   []float64    `json:"b"`
		A   []float64    `json:"a"`
		Sos [][6]float64 `json:"sos"`
	} `json:"output"`
}

// besselNorms is each of SciPy's norms as a dsp.BesselNorm
var besselNorms = map[string]dsp.BesselNorm{
	"mag":   dsp.BesselNormMagnitude,
	"phase": dsp.BesselNormPhase,
	"delay": dsp.BesselNormDelay,
}

// Bessel, measured on both architectures. The SciPy vectors are the exact
// design in each of SciPy's normalisations, computed with mpmath at 50 digits
// and rounded once, rather than scipy.signal.bessel's output, which is within
// 2.7e-14 of it. Both (b, a) and the sections are worst at 6.7e-15 on
// bandstop_order2_phase_wn0.01-0.99, whose wide edges are also Butter's and
// Cheby2's worst.
//
// MATLAB has no digital Bessel design: its vectors are besself at the
// prewarped edges through bilinear, for the phase normalisation, the only one
// besselap has, up to its highest order, 25. Worst at 1.6e-13 on the same
// bandstop's b, where MATLAB's own design is 1.7e-13 from the exact one.
// MATLAB's sos isn't compared, for the reason at butterMATLABSuite.
var (
	besselSciPySuite    = newSuite("Bessel", reference.SourceSciPy, "bessel", tolerance(1e-13), besselBAOutputs)
	besselSOSSciPySuite = newSuite("Bessel (SOS)", reference.SourceSciPy, "bessel", tolerance(1e-13), besselSOSOutputs)
	besselMATLABSuite   = newSuite("Bessel", reference.SourceMATLAB, "bessel", tolerance(1e-11), besselBAOutputs)
)

func designBessel(v besselVector) (f dsp.ZPK, err error) {
	wn, err := parseWn(v.Params.Wn)
	if err != nil {
		return f, err
	}
	band, err := bandFromVector(v.Params.BType, wn)
	if err != nil {
		return f, err
	}
	norm, ok := besselNorms[v.Params.Norm]
	if !ok {
		return f, fmt.Errorf("unknown norm %q", v.Params.Norm)
	}
	return dsp.Bessel(v.Params.Order, band, dsp.BesselNormalised(norm))
}

// besselBAOutputs pairs a design's (b, a) with a vector's
func besselBAOutputs(v besselVector) (outputs []Output, err error) {
	f, err := designBessel(v)
	if err != nil {
		return nil, err
	}
	b, a, err := f.BA()
	if err != nil {
		return nil, err
	}
	outputs = []Output{
		{Name: "b", Got: b, Want: v.Output.B},
		{Name: "a", Got: a, Want: v.Output.A},
	}
	return
}

// besselSOSOutputs pairs a design's sections with a vector's, compared as a
// set with the gain apart, as butterSOSSciPySuite explains
func besselSOSOutputs(v besselVector) (outputs []Output, err error) {
	f, err := designBessel(v)
	if err != nil {
		return nil, err
	}
	return matchedSOSOutputs(f, v.Output.Sos)
}
