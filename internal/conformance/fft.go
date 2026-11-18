package conformance

import (
	"github.com/cxd309/tukey/fft"
	"github.com/cxd309/tukey/internal/reference"
)

type fftVector struct {
	reference.Meta
	Params struct {
		X reference.Complex `json:"x"`
	} `json:"params"`
	Output struct {
		Y reference.Complex `json:"y"`
	} `json:"output"`
}

// forwardSciPySuite checks fft.Forward against scipy.fft.fft, at 23 lengths covering radix-2, mixed radix and Bluestein.
// Measured:
//   - radix-2 and mixed radix within 2 ULP
//   - Bluestein, at the primes 17, 97 and 1009, worst 4 ULP, 9.3e-16 at n = 1009
//
// Tolerance: about 10 times the worst.
var forwardSciPySuite = newSuite("fft.Forward", reference.SourceSciPy, "fft",
	reference.Tolerance{Rel: 1e-14, Scaled: 1e-14},
	forwardOutputs)

// forwardMATLABSuite checks fft.Forward against MATLAB's fft, which is FFTW.
// Measured: worst 4 ULP, 9.5e-16 at n = 1009, the same as against SciPy.
// Tolerance: the same as against SciPy.
var forwardMATLABSuite = newSuite("fft.Forward", reference.SourceMATLAB, "fft",
	reference.Tolerance{Rel: 1e-14, Scaled: 1e-14},
	forwardOutputs)

// inverseSciPySuite checks fft.Inverse against scipy.fft.ifft.
// Measured: worst 4 ULP, 8.5e-16 at n = 1009, the same as Forward.
// Tolerance: the same as Forward's.
var inverseSciPySuite = newSuite("fft.Inverse", reference.SourceSciPy, "ifft",
	reference.Tolerance{Rel: 1e-14, Scaled: 1e-14},
	inverseOutputs)

// inverseMATLABSuite checks fft.Inverse against MATLAB's ifft.
// Measured: worst 4 ULP, 8.3e-16 at n = 1009.
// Tolerance: the same as Forward's.
var inverseMATLABSuite = newSuite("fft.Inverse", reference.SourceMATLAB, "ifft",
	reference.Tolerance{Rel: 1e-14, Scaled: 1e-14},
	inverseOutputs)

// forwardOutputs pairs fft.Forward's output with a vector's.
func forwardOutputs(v fftVector) (outputs []Output, err error) {
	outputs, err = compareTransform(v, fft.Forward)
	return
}

// inverseOutputs pairs fft.Inverse's output with a vector's.
func inverseOutputs(v fftVector) (outputs []Output, err error) {
	outputs, err = compareTransform(v, fft.Inverse)
	return
}

// compareTransform runs a complex-to-complex transform on a vector's input, and pairs its output with the vector's.
func compareTransform(v fftVector, transform func([]complex128, ...fft.Option) ([]complex128, error)) (outputs []Output, err error) {
	x, err := v.Params.X.Values()
	if err != nil {
		return nil, err
	}
	want, err := v.Output.Y.Values()
	if err != nil {
		return nil, err
	}
	y, err := transform(x)
	if err != nil {
		return nil, err
	}
	outputs = []Output{{Name: "y", GotComplex: y, WantComplex: want}}
	return
}

type rfftVector struct {
	reference.Meta
	Params struct {
		X []float64 `json:"x"`
	} `json:"params"`
	Output struct {
		Y reference.Complex `json:"y"`
	} `json:"output"`
}

// forwardRealSciPySuite checks fft.ForwardReal against scipy.fft.rfft.
// Measured: worst 3 ULP, 5.8e-16 at n = 1009, by Bluestein,
// although SciPy uses a dedicated real-input algorithm where this uses the complex transform.
// Tolerance: the same as Forward's.
var forwardRealSciPySuite = newSuite("fft.ForwardReal", reference.SourceSciPy, "rfft",
	reference.Tolerance{Rel: 1e-14, Scaled: 1e-14},
	forwardRealOutputs)

// forwardRealMATLABSuite checks fft.ForwardReal against MATLAB's fft, cropped to its first n/2+1 bins,
// as MATLAB has no real-input transform.
// Measured: worst 3 ULP, 6.7e-16 at n = 1009.
// Tolerance: the same as Forward's.
var forwardRealMATLABSuite = newSuite("fft.ForwardReal", reference.SourceMATLAB, "rfft",
	reference.Tolerance{Rel: 1e-14, Scaled: 1e-14},
	forwardRealOutputs)

// forwardRealOutputs pairs fft.ForwardReal's output with a vector's.
func forwardRealOutputs(v rfftVector) (outputs []Output, err error) {
	want, err := v.Output.Y.Values()
	if err != nil {
		return nil, err
	}
	y, err := fft.ForwardReal(v.Params.X)
	if err != nil {
		return nil, err
	}
	outputs = []Output{{Name: "y", GotComplex: y, WantComplex: want}}
	return
}

type irfftVector struct {
	reference.Meta
	Params struct {
		X reference.Complex `json:"x"`
		N int               `json:"n"`
	} `json:"params"`
	Output struct {
		Y []float64 `json:"y"`
	} `json:"output"`
}

// inverseRealSciPySuite checks fft.InverseReal against scipy.fft.irfft.
// Measured: worst 2 ULP, 5.0e-16 at n = 1009.
// Tolerance: the same as Inverse's.
var inverseRealSciPySuite = newSuite("fft.InverseReal", reference.SourceSciPy, "irfft",
	reference.Tolerance{Rel: 1e-14, Scaled: 1e-14},
	inverseRealOutputs)

// inverseRealMATLABSuite checks fft.InverseReal against MATLAB's ifft(x, n, 'symmetric').
// Measured: worst 2 ULP, 5.3e-16 at n = 1009.
// Tolerance: the same as Inverse's.
var inverseRealMATLABSuite = newSuite("fft.InverseReal", reference.SourceMATLAB, "irfft",
	reference.Tolerance{Rel: 1e-14, Scaled: 1e-14},
	inverseRealOutputs)

// inverseRealOutputs pairs fft.InverseReal's output with a vector's.
func inverseRealOutputs(v irfftVector) (outputs []Output, err error) {
	x, err := v.Params.X.Values()
	if err != nil {
		return nil, err
	}
	y, err := fft.InverseReal(x, v.Params.N)
	if err != nil {
		return nil, err
	}
	outputs = []Output{{Name: "y", Got: y, Want: v.Output.Y}}
	return
}
