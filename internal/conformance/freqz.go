package conformance

import (
	"github.com/cxd309/tukey/dsp"
	"github.com/cxd309/tukey/internal/reference"
)

type freqzVector struct {
	reference.Meta
	Params struct {
		B []float64 `json:"b"`
		A []float64 `json:"a"`
		N int       `json:"n"`
	} `json:"params"`
	Output struct {
		W []float64         `json:"w"`
		H reference.Complex `json:"h"`
	} `json:"output"`
}

// freqzSuite checks dsp.Freqz against scipy.signal.freqz
// measured: worst 1.2e-15 (~6 ULP); Horner as SciPy, differing only in rounding of
// e^(-jw); fir_avg4, where SciPy uses an FFT instead, agrees as closely; set ~10x above
var freqzSuite = newSuite("Freqz", "freqz",
	reference.Tolerance{Rel: 1e-13, Scaled: 1e-13},
	func(v freqzVector) (outputs []Output, err error) {
		w, h, err := dsp.Freqz(v.Params.B, v.Params.A, v.Params.N)
		if err != nil {
			return nil, err
		}
		outputs, err = responseOutputs(w, h, v.Output.W, v.Output.H)
		return
	})

type sosFreqzVector struct {
	reference.Meta
	Params struct {
		SOS [][6]float64 `json:"sos"`
		N   int          `json:"n"`
	} `json:"params"`
	Output struct {
		W []float64         `json:"w"`
		H reference.Complex `json:"h"`
	} `json:"output"`
}

// sosFreqzSuite checks dsp.SOSFreqz against scipy.signal.freqz_sos
// measured: worst 3.7e-13 on butter8_lp0.01, whose poles sit so close to z = 1 that
// each section's denominator near DC is a near-cancellation, magnifying single-ULP
// differences in cos(w) between Go and NumPy; other filters within a few ULP;
// set ~10x above the worst
var sosFreqzSuite = newSuite("SOSFreqz", "freqz_sos",
	reference.Tolerance{Rel: 1e-11, Scaled: 1e-11},
	func(v sosFreqzVector) (outputs []Output, err error) {
		w, h, err := dsp.SOSFreqz(v.Params.SOS, v.Params.N)
		if err != nil {
			return nil, err
		}
		outputs, err = responseOutputs(w, h, v.Output.W, v.Output.H)
		return
	})

// responseOutputs pairs a frequency response with SciPy's: the frequencies,
// then the complex response
func responseOutputs(w []float64, h []complex128, wantW []float64, wantH reference.Complex) (outputs []Output, err error) {
	wantValues, err := wantH.Values()
	if err != nil {
		return nil, err
	}
	outputs = []Output{
		{Name: "w", Got: w, Want: wantW},
		{Name: "h", GotComplex: h, WantComplex: wantValues},
	}
	return
}
