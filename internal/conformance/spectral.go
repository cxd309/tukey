package conformance

import (
	"github.com/cxd309/tukey/dsp"
	"github.com/cxd309/tukey/internal/reference"
)

// complexRows is a complex matrix as stored in reference vector JSON: one row per segment
type complexRows struct {
	Re [][]float64 `json:"re"`
	Im [][]float64 `json:"im"`
}

// flattened is the matrix as one slice, row by row, so it's compared at one scale
func (c complexRows) flattened() (values []complex128) {
	for i, row := range c.Re {
		for j, re := range row {
			values = append(values, complex(re, c.Im[i][j]))
		}
	}
	return
}

type stftVector struct {
	reference.Meta
	Params struct {
		X        []float64 `json:"x"`
		FS       float64   `json:"fs"`
		Window   []float64 `json:"window"`
		Overlap  int       `json:"overlap"`
		NFFT     *int      `json:"nfft"` // null for the window's length
		TwoSided bool      `json:"two_sided"`
	} `json:"params"`
	Output struct {
		Frequencies []float64   `json:"frequencies"`
		Times       []float64   `json:"times"`
		Values      complexRows `json:"values"`
	} `json:"output"`
}

// stft measured: worst 4.2e-16 (2 ULP) against both, on tones_hamming255_nfft301,
// on both architectures. Set ~10x above.
var (
	stftSciPySuite  = newSuite("STFT", reference.SourceSciPy, "stft", tolerance(1e-14), stftOutputs)
	stftMATLABSuite = newSuite("STFT", reference.SourceMATLAB, "stft", tolerance(1e-14), stftOutputs)
)

// stftOutputs pairs dsp.STFT's frequencies, times and values with a vector's
func stftOutputs(v stftVector) (outputs []Output, err error) {
	p := v.Params
	opts := []dsp.SpectralOption{dsp.Overlap(p.Overlap)}
	if p.NFFT != nil {
		opts = append(opts, dsp.NFFT(*p.NFFT))
	}
	if p.TwoSided {
		opts = append(opts, dsp.TwoSided())
	}
	tf, err := dsp.STFT(p.X, p.FS, p.Window, opts...)
	if err != nil {
		return nil, err
	}
	frequencies := tf.Frequencies
	if p.TwoSided && v.Output.Frequencies[len(v.Output.Frequencies)-1] > p.FS/2 {
		// MATLAB labels the upper half above Nyquist, where tukey and SciPy label it negative
		frequencies = make([]float64, len(tf.Frequencies))
		for i, f := range tf.Frequencies {
			if f < 0 {
				f += p.FS
			}
			frequencies[i] = f
		}
	}
	outputs = []Output{
		{Name: "frequencies", Got: frequencies, Want: v.Output.Frequencies},
		{Name: "times", Got: tf.Times, Want: v.Output.Times},
		{Name: "values", GotComplex: flattenedValues(tf.Values), WantComplex: v.Output.Values.flattened()},
	}
	return
}

type stftFrequenciesVector struct {
	reference.Meta
	Params struct {
		X           []float64 `json:"x"`
		FS          float64   `json:"fs"`
		Window      []float64 `json:"window"`
		Overlap     int       `json:"overlap"`
		Frequencies []float64 `json:"frequencies"`
	} `json:"params"`
	Output struct {
		Frequencies []float64   `json:"frequencies"`
		Times       []float64   `json:"times"`
		Values      complexRows `json:"values"`
	} `json:"output"`
}

// stft_frequencies measured: worst 2.2e-16 (1 ULP) against SciPy, on tones_between_bins,
// and 2.4e-13 against MATLAB, on band_zoom, on both architectures. Set ~10x above.
// The SciPy side is summed exactly, with each twiddle's angle reduced exactly.
// MATLAB's frequency-vector path is that far from the exact sums itself, even at the
// FFT's own bins, where its FFT path is within 3e-16.
var (
	stftFrequenciesSciPySuite  = newSuite("STFT (AtFrequencies)", reference.SourceSciPy, "stft_frequencies", tolerance(1e-14), stftFrequenciesOutputs)
	stftFrequenciesMATLABSuite = newSuite("STFT (AtFrequencies)", reference.SourceMATLAB, "stft_frequencies", tolerance(1e-11), stftFrequenciesOutputs)
)

// stftFrequenciesOutputs pairs dsp.STFT with AtFrequencies against a vector's.
func stftFrequenciesOutputs(v stftFrequenciesVector) (outputs []Output, err error) {
	p := v.Params
	tf, err := dsp.STFT(p.X, p.FS, p.Window, dsp.Overlap(p.Overlap), dsp.AtFrequencies(p.Frequencies))
	if err != nil {
		return nil, err
	}
	outputs = []Output{
		{Name: "frequencies", Got: tf.Frequencies, Want: v.Output.Frequencies},
		{Name: "times", Got: tf.Times, Want: v.Output.Times},
		{Name: "values", GotComplex: flattenedValues(tf.Values), WantComplex: v.Output.Values.flattened()},
	}
	return
}

// flattenedValues is a TimeFrequency's values as one slice, segment by segment.
func flattenedValues(values [][]complex128) (flat []complex128) {
	for _, row := range values {
		flat = append(flat, row...)
	}
	return
}
