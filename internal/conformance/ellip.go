package conformance

import (
	"github.com/cxd309/tukey/dsp"
	"github.com/cxd309/tukey/internal/reference"
)

// Ellip, measured on both architectures. The SciPy vectors are the exact
// design, computed with mpmath at 50 digits and rounded once, rather than
// scipy.signal.ellip's output, which is within 1.5e-13 of it. Against them,
// (b, a) is worst at 3.9e-15 on lowpass_order10_rp1.0_rs40.0_wn0.3's b, and
// the sections at 3.0e-15. Against MATLAB it's worst at 9.9e-9 on
// lowpass_order4_rp0.01_rs80.0_wn0.3's b: MATLAB's ellipdeg forms 1 − k1² for
// the specification's small modulus k1 and loses its digits, where SciPy and
// tukey keep them. MATLAB's sos isn't compared, for the reason at
// butterMATLABSuite. The vectors reuse chebyVector, which has both dB params.
var (
	ellipSciPySuite    = newSuite("Ellip", reference.SourceSciPy, "ellip", tolerance(1e-13), chebyBAOutputs(designEllip))
	ellipSOSSciPySuite = newSuite("Ellip (SOS)", reference.SourceSciPy, "ellip", tolerance(1e-13), chebySOSOutputs(designEllip))
	ellipMATLABSuite   = newSuite("Ellip", reference.SourceMATLAB, "ellip", tolerance(1e-7), chebyBAOutputs(designEllip))
)

func designEllip(v chebyVector) (f dsp.ZPK, err error) {
	band, err := chebyBand(v)
	if err != nil {
		return f, err
	}
	return dsp.Ellip(v.Params.Order, v.Params.RippleDB, v.Params.AttenuationDB, band)
}
