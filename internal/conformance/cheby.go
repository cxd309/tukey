package conformance

import (
	"encoding/json"

	"github.com/cxd309/tukey/dsp"
	"github.com/cxd309/tukey/internal/reference"
)

type chebyVector struct {
	reference.Meta
	Params struct {
		Order         int             `json:"order"`
		RippleDB      float64         `json:"ripple_db"`      // cheby1 and ellip
		AttenuationDB float64         `json:"attenuation_db"` // cheby2 and ellip
		Wn            json.RawMessage `json:"wn"`
		BType         string          `json:"btype"`
	} `json:"params"`
	Output struct {
		B   []float64    `json:"b"`
		A   []float64    `json:"a"`
		Sos [][6]float64 `json:"sos"`
	} `json:"output"`
}

// Cheby1 and Cheby2, measured on both architectures. Against SciPy, Cheby1 is
// worst at 6.8e-15 on the 0.01 dB vector, where tukey computes the ripple
// factor by Expm1 and SciPy by 10^x − 1. Cheby2 is worst at 5.7e-14 on bandstop
// 0.01-0.99, also Butter's worst, whose wide edges make the poles
// ill-conditioned. Against MATLAB, Cheby1 is worst at 6.8e-13 on
// lowpass_order8_db0.5_wn0.01's tiny b, where MATLAB differs from SciPy by the
// same amount, and Cheby2 at 1.8e-14. MATLAB's sos isn't compared, for the
// reason at butterMATLABSuite.
var (
	cheby1SciPySuite    = newSuite("Cheby1", reference.SourceSciPy, "cheby1", tolerance(1e-13), chebyBAOutputs(designCheby1))
	cheby1SOSSciPySuite = newSuite("Cheby1 (SOS)", reference.SourceSciPy, "cheby1", tolerance(1e-13), chebySOSOutputs(designCheby1))
	cheby1MATLABSuite   = newSuite("Cheby1", reference.SourceMATLAB, "cheby1", tolerance(1e-11), chebyBAOutputs(designCheby1))
	cheby2SciPySuite    = newSuite("Cheby2", reference.SourceSciPy, "cheby2", tolerance(1e-12), chebyBAOutputs(designCheby2))
	cheby2SOSSciPySuite = newSuite("Cheby2 (SOS)", reference.SourceSciPy, "cheby2", tolerance(1e-12), chebySOSOutputs(designCheby2))
	cheby2MATLABSuite   = newSuite("Cheby2", reference.SourceMATLAB, "cheby2", tolerance(1e-12), chebyBAOutputs(designCheby2))
)

// chebyBand is the dsp.Band a vector describes.
func chebyBand(v chebyVector) (band dsp.Band, err error) {
	wn, err := parseWn(v.Params.Wn)
	if err != nil {
		return band, err
	}
	return bandFromVector(v.Params.BType, wn)
}

func designCheby1(v chebyVector) (f dsp.ZPK, err error) {
	band, err := chebyBand(v)
	if err != nil {
		return f, err
	}
	return dsp.Cheby1(v.Params.Order, v.Params.RippleDB, band)
}

func designCheby2(v chebyVector) (f dsp.ZPK, err error) {
	band, err := chebyBand(v)
	if err != nil {
		return f, err
	}
	return dsp.Cheby2(v.Params.Order, v.Params.AttenuationDB, band)
}

// chebyBAOutputs pairs a design's (b, a) with a vector's.
func chebyBAOutputs(design func(chebyVector) (dsp.ZPK, error)) func(chebyVector) ([]Output, error) {
	return func(v chebyVector) (outputs []Output, err error) {
		f, err := design(v)
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
}

// chebySOSOutputs pairs a design's sections with a vector's, compared as a
// set with the gain apart, as butterSOSSciPySuite explains.
func chebySOSOutputs(design func(chebyVector) (dsp.ZPK, error)) func(chebyVector) ([]Output, error) {
	return func(v chebyVector) (outputs []Output, err error) {
		f, err := design(v)
		if err != nil {
			return nil, err
		}
		return matchedSOSOutputs(f, v.Output.Sos)
	}
}
