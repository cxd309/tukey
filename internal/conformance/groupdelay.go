package conformance

import (
	"strings"

	"github.com/cxd309/tukey/dsp"
	"github.com/cxd309/tukey/internal/reference"
)

type groupDelayVector struct {
	reference.Meta
	Params struct {
		B []float64 `json:"b"`
		A []float64 `json:"a"`
		N int       `json:"n"`
	} `json:"params"`
	Output struct {
		W     []float64 `json:"w"`
		Delay []float64 `json:"delay"`
		// Singular lists the points where a zero or pole is on the unit
		// circle, whose delay is its limit from either side. The exact
		// vectors have it, MATLAB's don't.
		Singular []int `json:"singular"`
	} `json:"output"`
}

type sosGroupDelayVector struct {
	reference.Meta
	Params struct {
		Sos [][6]float64 `json:"sos"`
		N   int          `json:"n"`
	} `json:"params"`
	Output struct {
		W        []float64 `json:"w"`
		Delay    []float64 `json:"delay"`
		Singular []int     `json:"singular"`
	} `json:"output"`
}

// GroupDelay and SOSGroupDelay, measured on both architectures. The SciPy
// vectors are the exact group delay, computed with mpmath at 50 digits and
// rounded once, with a point on the unit circle's zero or pole at its limit
// from either side, rather than scipy.signal.group_delay's output, which is
// 3% off on butter8_lp0.01 and 0 or worse at those points. GroupDelay is
// worst at 7.2e-15 on ellip4_lp0.3, and SOSGroupDelay at 1.0e-14 on
// butter10_bp0.3-0.31, both at points kept in float64 by
// delayConditionLimit; butter8_lp0.01's (b, a), through big.Float, is exact.
//
// MATLAB's grpdelay loses digits near zeros and poles crowded against the
// unit circle, as SciPy's does. Its (b, a) vectors are split by how well
// conditioned they are: worst at 3.2e-11 on ellip4_lp0.3, and at 2.2e-5 on
// butter2_hp0.0008, near DC, among the two whose poles crowd against z = 1.
// butter8_lp0.01, which MATLAB has 7.3% off, isn't compared. Its sections
// are worst at 5.7e-9, on butter8_lp0.01's.
var (
	groupDelaySciPySuite    = newSuite("GroupDelay", reference.SourceSciPy, "group_delay", tolerance(1e-13), groupDelayOutputs)
	sosGroupDelaySciPySuite = newSuite("SOSGroupDelay", reference.SourceSciPy, "group_delay_sos", tolerance(1e-13), sosGroupDelayOutputs)

	groupDelayMATLABSuite = newSuiteOf("GroupDelay", reference.SourceMATLAB, "group_delay", tolerance(1e-9),
		groupDelayOutputs, func(file string) bool { return !illConditionedDelay(file) && !matlabDelayUncompared(file) })
	groupDelayIllConditionedMATLABSuite = newSuiteOf("GroupDelay (ill-conditioned)", reference.SourceMATLAB, "group_delay", tolerance(1e-3),
		groupDelayOutputs, illConditionedDelay)
	sosGroupDelayMATLABSuite = newSuite("SOSGroupDelay", reference.SourceMATLAB, "group_delay_sos", tolerance(1e-7), sosGroupDelayOutputs)
)

// illConditionedDelay accepts the (b, a) vectors whose poles crowd against
// z = 1, where MATLAB's grpdelay is only good to a few digits
func illConditionedDelay(file string) bool {
	return strings.HasPrefix(file, "butter2_lp0.001_") || strings.HasPrefix(file, "butter2_hp0.0008_")
}

// matlabDelayUncompared accepts the vectors MATLAB's grpdelay is too far off
// to compare
func matlabDelayUncompared(file string) bool {
	return strings.HasPrefix(file, "butter8_lp0.01_")
}

// groupDelayOutputs pairs dsp.GroupDelay's frequencies and delay with a vector's
func groupDelayOutputs(v groupDelayVector) (outputs []Output, err error) {
	w, delay, err := dsp.GroupDelay(v.Params.B, v.Params.A, v.Params.N)
	if err != nil {
		return nil, err
	}
	outputs = []Output{
		{Name: "w", Got: w, Want: v.Output.W},
		{Name: "delay", Got: delay, Want: v.Output.Delay},
	}
	return
}

// sosGroupDelayOutputs pairs dsp.SOSGroupDelay's frequencies and delay with a vector's
func sosGroupDelayOutputs(v sosGroupDelayVector) (outputs []Output, err error) {
	sos := make(dsp.SOS, len(v.Params.Sos))
	for i, section := range v.Params.Sos {
		sos[i] = section
	}
	w, delay, err := dsp.SOSGroupDelay(sos, v.Params.N)
	if err != nil {
		return nil, err
	}
	outputs = []Output{
		{Name: "w", Got: w, Want: v.Output.W},
		{Name: "delay", Got: delay, Want: v.Output.Delay},
	}
	return
}
