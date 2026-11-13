package conformance

import (
	"fmt"
	"math/cmplx"

	"github.com/cxd309/tukey/dsp"
	"github.com/cxd309/tukey/internal/reference"
)

type notchVector struct {
	reference.Meta
	Params struct {
		W0 float64 `json:"w0"`
		Q  float64 `json:"q"`
	} `json:"params"`
	Output struct {
		B []float64         `json:"b"`
		A []float64         `json:"a"`
		Z reference.Complex `json:"z"`
		P reference.Complex `json:"p"`
		K float64           `json:"k"`
	} `json:"output"`
}

// IIRNotch and IIRPeak, measured on both architectures. The SciPy vectors
// are the exact design, computed with mpmath at 50 digits and rounded once,
// rather than scipy.signal's output: its notch agrees to 3.8e-16, but its
// peak's gain is 1 − 1/(1 + β), which cancels for a narrow bandwidth, off by
// 6.1e-13 at w0 = 0.002, Q = 30. Both are worst at 4.7e-16 (2 ULP), on
// w00.3_q1000.0's a, whose poles are hard against the unit circle.
//
// Not yet checked against MATLAB: its iirnotch and iirpeak are in the DSP
// System Toolbox.
var (
	iirNotchSciPySuite = newSuite("IIRNotch", reference.SourceSciPy, "iirnotch", tolerance(1e-14), notchOutputs(dsp.IIRNotch))
	iirPeakSciPySuite  = newSuite("IIRPeak", reference.SourceSciPy, "iirpeak", tolerance(1e-14), notchOutputs(dsp.IIRPeak))
)

// notchOutputs pairs a design's (b, a), zeros, poles and gain with a vector's.
// The zeros and poles are compared as sets, each wanted root against the
// nearest designed one.
func notchOutputs(design func(w0, q float64) (dsp.ZPK, error)) func(notchVector) ([]Output, error) {
	return func(v notchVector) (outputs []Output, err error) {
		f, err := design(v.Params.W0, v.Params.Q)
		if err != nil {
			return nil, err
		}
		b, a, err := f.BA()
		if err != nil {
			return nil, err
		}
		zeros, err := v.Output.Z.Values()
		if err != nil {
			return nil, err
		}
		poles, err := v.Output.P.Values()
		if err != nil {
			return nil, err
		}
		gotZeros, err := nearestRoots(f.Zeros, zeros)
		if err != nil {
			return nil, fmt.Errorf("zeros: %w", err)
		}
		gotPoles, err := nearestRoots(f.Poles, poles)
		if err != nil {
			return nil, fmt.Errorf("poles: %w", err)
		}
		outputs = []Output{
			{Name: "b", Got: b, Want: v.Output.B},
			{Name: "a", Got: a, Want: v.Output.A},
			{Name: "z", GotComplex: gotZeros, WantComplex: zeros},
			{Name: "p", GotComplex: gotPoles, WantComplex: poles},
			{Name: "k", Got: []float64{f.Gain}, Want: []float64{v.Output.K}},
		}
		return
	}
}

// nearestRoots is got reordered to line up with want: each wanted root takes
// the nearest of the got roots not yet taken
func nearestRoots(got, want []complex128) (ordered []complex128, err error) {
	if len(got) != len(want) {
		return nil, fmt.Errorf("got %d roots, want %d", len(got), len(want))
	}
	taken := make([]bool, len(got))
	ordered = make([]complex128, len(want))
	for i, w := range want {
		best := -1
		for j, g := range got {
			if !taken[j] && (best < 0 || cmplx.Abs(g-w) < cmplx.Abs(got[best]-w)) {
				best = j
			}
		}
		taken[best] = true
		ordered[i] = got[best]
	}
	return
}
