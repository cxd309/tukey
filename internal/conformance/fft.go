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
func compareTransform(v fftVector, transform func([]complex128) ([]complex128, error)) (outputs []Output, err error) {
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
