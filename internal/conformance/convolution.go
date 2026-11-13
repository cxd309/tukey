package conformance

import (
	"github.com/cxd309/tukey/dsp"
	"github.com/cxd309/tukey/internal/reference"
)

type convolveVector struct {
	reference.Meta
	Params struct {
		A []float64 `json:"a"`
		V []float64 `json:"v"`
	} `json:"params"`
	Output struct {
		Full  []float64 `json:"full"`
		Same  []float64 `json:"same"`
		Valid []float64 `json:"valid"`
	} `json:"output"`
}

type correlateVector struct {
	reference.Meta
	Params struct {
		A []float64 `json:"a"`
		V []float64 `json:"v"`
	} `json:"params"`
	Output struct {
		Full      []float64 `json:"full"`
		Same      []float64 `json:"same"`  // SciPy only
		Valid     []float64 `json:"valid"` // SciPy only
		FullLags  []int     `json:"full_lags"`
		SameLags  []int     `json:"same_lags"`  // SciPy only
		ValidLags []int     `json:"valid_lags"` // SciPy only
	} `json:"output"`
}

// Convolve and Correlate, measured on both architectures. The SciPy side is
// exact, each output its exact sum of exact products, so both methods are
// judged against the true answer: direct matches it exactly, and FFT is
// within 5.5e-16 of the largest output (fir_lowpass). MATLAB's conv is
// direct, and itself up to 1.6e-15 off the exact answer on long_both, where
// FFT is 1.8e-15 off it. Its xcorr is by FFT, 2.4e-4 off the exact answer on
// mixed_magnitudes' smallest outputs, though within 6.3e-16 of the largest,
// the scale every comparison here is measured against.
var (
	convolveDirectSciPySuite   = newSuite("Convolve (direct)", reference.SourceSciPy, "convolve", tolerance(1e-14), convolveOutputs(dsp.ConvolveDirect, reference.SourceSciPy))
	convolveDirectMATLABSuite  = newSuite("Convolve (direct)", reference.SourceMATLAB, "convolve", tolerance(1e-13), convolveOutputs(dsp.ConvolveDirect, reference.SourceMATLAB))
	convolveFFTSciPySuite      = newSuite("Convolve (FFT)", reference.SourceSciPy, "convolve", tolerance(1e-14), convolveOutputs(dsp.ConvolveFFT, reference.SourceSciPy))
	convolveFFTMATLABSuite     = newSuite("Convolve (FFT)", reference.SourceMATLAB, "convolve", tolerance(1e-13), convolveOutputs(dsp.ConvolveFFT, reference.SourceMATLAB))
	correlateDirectSciPySuite  = newSuite("Correlate (direct)", reference.SourceSciPy, "correlate", tolerance(1e-14), correlateOutputs(dsp.ConvolveDirect, reference.SourceSciPy))
	correlateDirectMATLABSuite = newSuite("Correlate (direct)", reference.SourceMATLAB, "correlate", tolerance(1e-14), correlateOutputs(dsp.ConvolveDirect, reference.SourceMATLAB))
	correlateFFTSciPySuite     = newSuite("Correlate (FFT)", reference.SourceSciPy, "correlate", tolerance(1e-14), correlateOutputs(dsp.ConvolveFFT, reference.SourceSciPy))
	correlateFFTMATLABSuite    = newSuite("Correlate (FFT)", reference.SourceMATLAB, "correlate", tolerance(1e-14), correlateOutputs(dsp.ConvolveFFT, reference.SourceMATLAB))
)

// convolveMode is one of a vector's outputs: its name, the option selecting
// it, and the expected values.
type convolveMode struct {
	name string
	opts []dsp.ConvolveOption
	want []float64
}

// comparedModes is the modes a vector of na and nv samples is compared in.
// A mode the vector doesn't have (MATLAB's xcorr has no "same" or "valid")
// is skipped. So are MATLAB's "same" for an even nv, which starts one output
// later than tukey's, and its "valid" for nv > na, which is empty where
// tukey's swaps the inputs, both documented divergences.
func comparedModes(source reference.Source, na, nv int, full, same, valid []float64) (modes []convolveMode) {
	matlab := source == reference.SourceMATLAB
	candidates := []struct {
		mode     convolveMode
		compared bool
	}{
		{convolveMode{"full", nil, full}, true},
		{convolveMode{"same", []dsp.ConvolveOption{dsp.ConvolveSame()}, same}, !matlab || nv%2 == 1},
		{convolveMode{"valid", []dsp.ConvolveOption{dsp.ConvolveValid()}, valid}, !matlab || na >= nv},
	}
	for _, c := range candidates {
		if c.compared && c.mode.want != nil {
			modes = append(modes, c.mode)
		}
	}
	return
}

// convolveOutputs pairs dsp.Convolve by method, in each compared mode, with a vector's.
func convolveOutputs(method dsp.ConvolveMethod, source reference.Source) func(v convolveVector) ([]Output, error) {
	return func(v convolveVector) (outputs []Output, err error) {
		a, k := v.Params.A, v.Params.V
		for _, m := range comparedModes(source, len(a), len(k), v.Output.Full, v.Output.Same, v.Output.Valid) {
			c, err := dsp.Convolve(a, k, append(m.opts, dsp.ConvolveBy(method))...)
			if err != nil {
				return nil, err
			}
			outputs = append(outputs, Output{Name: m.name, Got: c, Want: m.want})
		}
		return
	}
}

// correlateOutputs pairs dsp.Correlate by method, in each compared mode, with
// a vector's values and lags. Lags are compared as floats, exactly.
func correlateOutputs(method dsp.ConvolveMethod, source reference.Source) func(v correlateVector) ([]Output, error) {
	return func(v correlateVector) (outputs []Output, err error) {
		a, k := v.Params.A, v.Params.V
		lags := map[string][]int{"full": v.Output.FullLags, "same": v.Output.SameLags, "valid": v.Output.ValidLags}
		for _, m := range comparedModes(source, len(a), len(k), v.Output.Full, v.Output.Same, v.Output.Valid) {
			c, err := dsp.Correlate(a, k, append(m.opts, dsp.ConvolveBy(method))...)
			if err != nil {
				return nil, err
			}
			outputs = append(outputs,
				Output{Name: m.name, Got: c.Values, Want: m.want},
				Output{Name: m.name + "_lags", Got: asFloats(c.Lags), Want: asFloats(lags[m.name])},
			)
		}
		return
	}
}

// asFloats is ints as float64s, to compare with Output.
func asFloats(ints []int) (floats []float64) {
	floats = make([]float64, len(ints))
	for i, n := range ints {
		floats[i] = float64(n)
	}
	return
}
