package conformance

import (
	"fmt"
	"slices"

	"github.com/cxd309/tukey/dsp"
	"github.com/cxd309/tukey/internal/reference"
)

// complexRows is a complex matrix as stored in reference vector JSON, one row per segment.
type complexRows struct {
	Re [][]float64 `json:"re"`
	Im [][]float64 `json:"im"`
}

// flattened is the matrix as one slice, row by row, so it's compared at one scale.
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
		NFFT     *int      `json:"nfft"` // or null for the window's length
		TwoSided bool      `json:"two_sided"`
	} `json:"params"`
	Output struct {
		Frequencies []float64   `json:"frequencies"`
		Times       []float64   `json:"times"`
		Values      complexRows `json:"values"`
	} `json:"output"`
}

// The stft suites.
// Measured: worst 2 ULP, 4.2e-16, against both, on tones_hamming255_nfft301, on both architectures.
// Tolerance: about 10 times the worst.
var (
	stftSciPySuite  = newSuite("STFT", reference.SourceSciPy, "stft", tolerance(1e-14), stftOutputs)
	stftMATLABSuite = newSuite("STFT", reference.SourceMATLAB, "stft", tolerance(1e-14), stftOutputs)
)

// stftOutputs pairs dsp.STFT's frequencies, times and values with a vector's.
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

// The stft_frequencies suites.
// Measured: worst 1 ULP, 2.2e-16, against SciPy, on tones_between_bins,
// and 2.4e-13 against MATLAB, on band_zoom, on both architectures.
// Tolerance: about 10 times the worst.
// The SciPy side is summed exactly, with each twiddle's angle reduced exactly.
// MATLAB's frequency-vector path is that far from the exact sums itself, even at the FFT's own bins,
// where its FFT path is within 3e-16.
var (
	stftFrequenciesSciPySuite  = newSuite("STFT (AtFrequencies)", reference.SourceSciPy, "stft_frequencies", tolerance(1e-14), stftFrequenciesOutputs)
	stftFrequenciesMATLABSuite = newSuite("STFT (AtFrequencies)", reference.SourceMATLAB, "stft_frequencies", tolerance(1e-11), stftFrequenciesOutputs)
)

// stftFrequenciesOutputs pairs dsp.STFT, with AtFrequencies, with a vector's.
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

type spectrogramVector struct {
	reference.Meta
	Params struct {
		X       []float64 `json:"x"`
		FS      float64   `json:"fs"`
		Window  []float64 `json:"window"`
		Overlap int       `json:"overlap"`
		NFFT    *int      `json:"nfft"`
		Scaling string    `json:"scaling"`
	} `json:"params"`
	Output struct {
		Frequencies []float64   `json:"frequencies"`
		Times       []float64   `json:"times"`
		Values      [][]float64 `json:"values"`
	} `json:"output"`
}

// The spectrogram suites.
// Measured: worst 3 ULP, 7.3e-16, against both, on chirp_hamming255_nfft301, on both architectures.
// Tolerance: about 10 times the worst.
var (
	spectrogramSciPySuite  = newSuite("Spectrogram", reference.SourceSciPy, "spectrogram", tolerance(1e-14), spectrogramOutputs)
	spectrogramMATLABSuite = newSuite("Spectrogram", reference.SourceMATLAB, "spectrogram", tolerance(1e-14), spectrogramOutputs)
)

// spectrogramOutputs pairs dsp.Spectrogram's frequencies, times and power with a vector's.
func spectrogramOutputs(v spectrogramVector) (outputs []Output, err error) {
	p := v.Params
	opts := append(powerOptions(p.Scaling, p.NFFT), dsp.Overlap(p.Overlap))
	tf, err := dsp.Spectrogram(p.X, p.FS, p.Window, opts...)
	if err != nil {
		return nil, err
	}
	outputs = []Output{
		{Name: "frequencies", Got: tf.Frequencies, Want: v.Output.Frequencies},
		{Name: "times", Got: tf.Times, Want: v.Output.Times},
		{Name: "values", Got: slices.Concat(tf.Values...), Want: slices.Concat(v.Output.Values...)},
	}
	return
}

type welchVector struct {
	reference.Meta
	Params struct {
		X       []float64 `json:"x"`
		FS      float64   `json:"fs"`
		Window  []float64 `json:"window"`
		Overlap int       `json:"overlap"`
		NFFT    *int      `json:"nfft"`
		Scaling string    `json:"scaling"`
		Average string    `json:"average"`
		Detrend string    `json:"detrend"`
	} `json:"params"`
	Output struct {
		Frequencies []float64 `json:"frequencies"`
		Power       []float64 `json:"power"`
	} `json:"output"`
}

// The welch suites.
// Measured: worst 5 ULP, 1.1e-15, against SciPy, and 6 ULP, 1.3e-15, against MATLAB on arm64,
// on chirp_hamming255_nfft512, on both architectures.
// Tolerance: about 10 times the worst.
// The median average, with its bias correction, and per-segment detrending match SciPy too.
// MATLAB's pwelch has neither, so those cases are SciPy only.
var (
	welchSciPySuite  = newSuite("Welch", reference.SourceSciPy, "welch", tolerance(1e-13), welchOutputs)
	welchMATLABSuite = newSuite("Welch", reference.SourceMATLAB, "welch", tolerance(1e-13), welchOutputs)
)

// welchOutputs pairs dsp.Welch's frequencies and power with a vector's.
func welchOutputs(v welchVector) (outputs []Output, err error) {
	p := v.Params
	opts := append(powerOptions(p.Scaling, p.NFFT), dsp.Overlap(p.Overlap))
	switch p.Average {
	case "mean":
	case "median":
		opts = append(opts, dsp.AverageMedian())
	default:
		return nil, fmt.Errorf("unknown average %q", p.Average)
	}
	switch p.Detrend {
	case "none":
	case "constant":
		opts = append(opts, dsp.SegmentDetrend(dsp.DetrendConstant()))
	case "linear":
		opts = append(opts, dsp.SegmentDetrend(dsp.DetrendLinear()))
	default:
		return nil, fmt.Errorf("unknown detrend %q", p.Detrend)
	}
	spectrum, err := dsp.Welch(p.X, p.FS, p.Window, opts...)
	if err != nil {
		return nil, err
	}
	outputs = []Output{
		{Name: "frequencies", Got: spectrum.Frequencies, Want: v.Output.Frequencies},
		{Name: "power", Got: spectrum.Power, Want: v.Output.Power},
	}
	return
}

type periodogramVector struct {
	reference.Meta
	Params struct {
		X       []float64 `json:"x"`
		FS      float64   `json:"fs"`
		Window  []float64 `json:"window"`
		NFFT    *int      `json:"nfft"`
		Scaling string    `json:"scaling"`
	} `json:"params"`
	Output struct {
		Frequencies []float64 `json:"frequencies"`
		Power       []float64 `json:"power"`
	} `json:"output"`
}

// The periodogram suites.
// Measured: worst 5 ULP, 1.1e-15, against SciPy, on tones_hann_spectrum,
// and 2 ULP, 5.0e-16, against MATLAB, on tones_boxcar, on both architectures.
// Tolerance: about 10 times the worst.
var (
	periodogramSciPySuite  = newSuite("Periodogram", reference.SourceSciPy, "periodogram", tolerance(1e-13), periodogramOutputs)
	periodogramMATLABSuite = newSuite("Periodogram", reference.SourceMATLAB, "periodogram", tolerance(1e-14), periodogramOutputs)
)

// periodogramOutputs pairs dsp.Periodogram's frequencies and power with a vector's.
func periodogramOutputs(v periodogramVector) (outputs []Output, err error) {
	p := v.Params
	spectrum, err := dsp.Periodogram(p.X, p.FS, p.Window, powerOptions(p.Scaling, p.NFFT)...)
	if err != nil {
		return nil, err
	}
	outputs = []Output{
		{Name: "frequencies", Got: spectrum.Frequencies, Want: v.Output.Frequencies},
		{Name: "power", Got: spectrum.Power, Want: v.Output.Power},
	}
	return
}

// powerOptions is a vector's scaling and NFFT as SpectralOptions.
// An unknown scaling is left as the default, density.
func powerOptions(scaling string, nfft *int) (opts []dsp.SpectralOption) {
	if scaling == "spectrum" {
		opts = append(opts, dsp.ScalingSpectrum())
	}
	if nfft != nil {
		opts = append(opts, dsp.NFFT(*nfft))
	}
	return
}
