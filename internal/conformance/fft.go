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

// forwardSuite checks fft.Forward against scipy.fft.fft at 23 lengths covering
// radix-2, mixed radix and Bluestein
// measured: radix-2 and mixed radix within 2 ULP; Bluestein (primes 17, 97, 1009)
// worst 4 ULP (9.3e-16 at n = 1009); set ~10x above the worst
var forwardSuite = newSuite("fft.Forward", "fft",
	reference.Tolerance{Rel: 1e-14, Scaled: 1e-14},
	func(v fftVector) (outputs []Output, err error) {
		outputs, err = compareTransform(v, fft.Forward)
		return
	})

// inverseSuite checks fft.Inverse against scipy.fft.ifft
// measured: as Forward, worst 4 ULP (8.5e-16 at n = 1009); set to match
var inverseSuite = newSuite("fft.Inverse", "ifft",
	reference.Tolerance{Rel: 1e-14, Scaled: 1e-14},
	func(v fftVector) (outputs []Output, err error) {
		outputs, err = compareTransform(v, fft.Inverse)
		return
	})

// compareTransform runs a complex-to-complex transform on a vector's input
// and pairs its output with SciPy's
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

// forwardRealSuite checks fft.ForwardReal against scipy.fft.rfft
// measured: worst 3 ULP (5.8e-16 at n = 1009, Bluestein) despite SciPy using a
// dedicated real-input algorithm where this uses the complex transform; set to
// match Forward
var forwardRealSuite = newSuite("fft.ForwardReal", "rfft",
	reference.Tolerance{Rel: 1e-14, Scaled: 1e-14},
	func(v rfftVector) (outputs []Output, err error) {
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
	})

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

// inverseRealSuite checks fft.InverseReal against scipy.fft.irfft
// measured: worst 2 ULP (5.0e-16 at n = 1009); set to match Inverse
var inverseRealSuite = newSuite("fft.InverseReal", "irfft",
	reference.Tolerance{Rel: 1e-14, Scaled: 1e-14},
	func(v irfftVector) (outputs []Output, err error) {
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
	})
