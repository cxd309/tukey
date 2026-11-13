package conformance

import (
	"encoding/json"

	"github.com/cxd309/tukey/dsp"
	"github.com/cxd309/tukey/internal/reference"
)

type firWinVector struct {
	reference.Meta
	Params struct {
		Window []float64       `json:"window"`
		Wn     json.RawMessage `json:"wn"`
		BType  string          `json:"btype"`
		Scale  bool            `json:"scale"`
	} `json:"params"`
	Output struct {
		H []float64 `json:"h"`
	} `json:"output"`
}

// FIRWin, measured on both architectures. The SciPy vectors are the exact
// design from the window's samples, computed with mpmath at 50 digits and
// rounded once, rather than scipy.signal.firwin's output, which is within
// 3.0e-15 of it. Worst at 3.5e-16 (2 ULP) on lowpass_hamming_n501_wn0.01,
// whose outer taps are furthest from the centre. MATLAB's are fir1 with the
// same window, skipping the one-tap case, which fir1 can't design: worst at
// 3.2e-15 on the narrow bandpass_hamming_n51_wn0.3-0.31, where fir1 is
// itself 3.2e-15 from the exact design.
var (
	firWinSciPySuite  = newSuite("FIRWin", reference.SourceSciPy, "firwin", tolerance(1e-14), firWinOutputs)
	firWinMATLABSuite = newSuite("FIRWin", reference.SourceMATLAB, "firwin", tolerance(1e-13), firWinOutputs)
)

// firWinOutputs pairs dsp.FIRWin's taps with a vector's
func firWinOutputs(v firWinVector) (outputs []Output, err error) {
	wn, err := parseWn(v.Params.Wn)
	if err != nil {
		return nil, err
	}
	band, err := bandFromVector(v.Params.BType, wn)
	if err != nil {
		return nil, err
	}
	var opts []dsp.FIROption
	if !v.Params.Scale {
		opts = append(opts, dsp.NoScale())
	}
	h, err := dsp.FIRWin(v.Params.Window, band, opts...)
	if err != nil {
		return nil, err
	}
	outputs = []Output{{Name: "h", Got: h, Want: v.Output.H}}
	return
}
