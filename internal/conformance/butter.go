package conformance

import (
	"encoding/json"
	"fmt"

	"github.com/cxd309/tukey/dsp"
	"github.com/cxd309/tukey/internal/reference"
)

type butterVector struct {
	reference.Meta
	Params struct {
		Order int             `json:"order"`
		Wn    json.RawMessage `json:"wn"`
		BType string          `json:"btype"`
	} `json:"params"`
	Output struct {
		B   []float64    `json:"b"`
		A   []float64    `json:"a"`
		Sos [][6]float64 `json:"sos"`
	} `json:"output"`
}

// butterSciPySuite checks dsp.Butter against scipy.signal.butter.
// Measured: most cases within a few ULP, about 1e-16, and worst about 7e-15 on bandstop 0.01-0.99,
// where prewarping near Nyquist makes the pole quadratic ill-conditioned.
// Tolerance: about 10 times the worst, for headroom across platforms.
var butterSciPySuite = newSuite("Butter", reference.SourceSciPy, "butter",
	reference.Tolerance{Rel: 1e-13, Scaled: 1e-13},
	butterOutputs)

// butterMATLABSuite checks dsp.Butter's (b, a) against MATLAB's butter.
// Its sos isn't compared, as zp2sos factorises odd orders and equal-radius pole pairs differently from SciPy's zpk2sos,
// giving the same filter in different sections.
// Measured: most cases within a few ULP, and worst 533 ULP, 1.2e-13, on lowpass_order8_wn0.01's b, of about 1e-14.
// There MATLAB's own b differs from SciPy's by the same amount, while tukey matches SciPy to 7e-15,
// so the difference is in MATLAB's design, not tukey's.
// Tolerance: about 10 times the worst.
var butterMATLABSuite = newSuite("Butter", reference.SourceMATLAB, "butter",
	reference.Tolerance{Rel: 1e-11, Scaled: 1e-11},
	butterOutputs)

// butterOutputs pairs dsp.Butter's (b, a) with a vector's.
func butterOutputs(v butterVector) (outputs []Output, err error) {
	wn, err := parseWn(v.Params.Wn)
	if err != nil {
		return nil, err
	}
	band, err := bandFromVector(v.Params.BType, wn)
	if err != nil {
		return nil, err
	}
	f, err := dsp.Butter(v.Params.Order, band)
	if err != nil {
		return nil, err
	}
	b, a, err := f.BA()
	if err != nil {
		return nil, err
	}
	// no Scale, as coefficients are their own scale
	outputs = []Output{
		{Name: "b", Got: b, Want: v.Output.B},
		{Name: "a", Got: a, Want: v.Output.A},
	}
	return
}

// parseWn handles the vectors' wn field,
// which is a scalar for a lowpass or highpass, and a [low, high] pair for a bandpass or bandstop.
func parseWn(raw json.RawMessage) (wn []float64, err error) {
	var scalar float64
	if err := json.Unmarshal(raw, &scalar); err == nil {
		return []float64{scalar}, nil
	}
	var pair []float64
	if err := json.Unmarshal(raw, &pair); err == nil {
		return pair, nil
	}
	return nil, fmt.Errorf("wn is neither a float nor a [low, high] pair: %s", raw)
}

// bandFromVector builds the dsp.Band a vector describes.
func bandFromVector(btype string, wn []float64) (band dsp.Band, err error) {
	switch btype {
	case "lowpass":
		band = dsp.Lowpass(wn[0])
	case "highpass":
		band = dsp.Highpass(wn[0])
	case "bandpass":
		band = dsp.Bandpass(wn[0], wn[1])
	case "bandstop":
		band = dsp.Bandstop(wn[0], wn[1])
	default:
		err = fmt.Errorf("unknown btype %q", btype)
	}
	return
}
