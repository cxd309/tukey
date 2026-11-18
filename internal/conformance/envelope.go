package conformance

import (
	"fmt"

	"github.com/cxd309/tukey/dsp"
	"github.com/cxd309/tukey/internal/reference"
)

type envelopeVector struct {
	reference.Meta
	Params struct {
		X    []float64 `json:"x"`
		Mode string    `json:"mode"`
		N    *int      `json:"n"` // the RMS window, or null for analytic
	} `json:"params"`
	Output struct {
		Upper []float64 `json:"upper"`
		Lower []float64 `json:"lower"`
	} `json:"output"`
}

// The envelope suites.
// Measured:
//   - worst 5 ULP, 1.0e-15, against SciPy, on sines_analytic on amd64, and 7.4e-16 on arm64
//   - worst 20 ULP, 4.5e-15, against MATLAB, on impulse_rms250, on both architectures
//
// Tolerance: about 10 times the worst.
//
// The compensated mean took sines_analytic from 3 to 5 ULP on amd64, although it's the exactly rounded mean,
// where the plain sum's was 4.2e-17 off.
// That small change in the centred input moves where the FFT's roundings land,
// so the difference from NumPy is rounding noise, not lost accuracy.
// The SciPy side is MATLAB's definition built from SciPy and NumPy, which MATLAB's own vectors confirm.
// The gap to MATLAB is MATLAB's own rounding.
// Where one large square sits among many small ones, MATLAB's RMS envelope loses the small ones' low bits,
// the same as a plain sum, so on impulse_rms250 it's 4.5e-15 from the exactly summed answer, by math.fsum,
// where the SciPy side, summed pairwise by NumPy, is 6.6e-16.
// dsp.Envelope's compensated sum matches the SciPy side.
var (
	envelopeSciPySuite  = newSuite("Envelope", reference.SourceSciPy, "envelope", tolerance(1e-14), envelopeOutputs)
	envelopeMATLABSuite = newSuite("Envelope", reference.SourceMATLAB, "envelope", tolerance(1e-13), envelopeOutputs)
)

// envelopeOutputs pairs dsp.Envelope's upper and lower envelopes with a vector's.
func envelopeOutputs(v envelopeVector) (outputs []Output, err error) {
	var method dsp.EnvelopeMethod
	switch v.Params.Mode {
	case "analytic":
		method = dsp.EnvelopeAnalytic()
	case "rms":
		if v.Params.N == nil {
			return nil, fmt.Errorf("rms envelope without a window")
		}
		method = dsp.EnvelopeRMS(*v.Params.N)
	default:
		return nil, fmt.Errorf("unknown envelope mode %q", v.Params.Mode)
	}
	upper, lower, err := dsp.Envelope(v.Params.X, method)
	if err != nil {
		return nil, err
	}
	outputs = []Output{
		{Name: "upper", Got: upper, Want: v.Output.Upper},
		{Name: "lower", Got: lower, Want: v.Output.Lower},
	}
	return
}
