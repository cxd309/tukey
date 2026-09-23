package conformance

import (
	"encoding/json"
	"fmt"

	"github.com/cxd309/godsp/dsp"
	"github.com/cxd309/godsp/internal/reference"
)

type butterVector struct {
	reference.Meta
	Params struct {
		Order int             `json:"order"`
		Wn    json.RawMessage `json:"wn"`
		BType string          `json:"btype"`
	} `json:"params"`
	Output struct {
		B []float64 `json:"b"`
		A []float64 `json:"a"`
	} `json:"output"`
}

// butterSuite checks dsp.Butter against scipy.signal.butter
// measured: most cases agree within a few ULPs (~1e-16)
// worst is bandstop 0.01-0.99 at ~7e-15
// where prewarp near Nyquist makes the pole quadratic ill-conditioned;
// set ~10x above the worst for headroom across platforms
var butterSuite = newSuite("Butter", "butter",
	reference.Tolerance{Rel: 1e-13, Abs: 1e-13},
	func(v butterVector) (outputs []Output, err error) {
		wn, err := parseWn(v.Params.Wn)
		if err != nil {
			return nil, err
		}
		band, err := bandFromVector(v.Params.BType, wn)
		if err != nil {
			return nil, err
		}
		b, a, err := dsp.Butter(v.Params.Order, band)
		if err != nil {
			return nil, err
		}
		outputs = []Output{{"b", b, v.Output.B}, {"a", a, v.Output.A}}
		return
	})

// parseWn handles the vectors' wn field, which is a scalar for
// lowpass/highpass and a [low, high] pair for bandpass/bandstop
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

// bandFromVector builds the dsp.Band a vector describes
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
