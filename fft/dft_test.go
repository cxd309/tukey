package fft

import (
	"testing"

	"github.com/cxd309/tukey/internal/reference"
	"github.com/stretchr/testify/require"
)

// dft is the discrete Fourier transform straight from its definition, in O(n²):
//
//	y[k] = Σ x[j]·e^(-2πi·jk/n)
//
// It's far too slow for real use, and is the reference the fast algorithms are checked against.
func dft(x []complex128) (y []complex128) {
	n := len(x)
	y = make([]complex128, n)
	for k := range y {
		var sum complex128
		for j, xj := range x {
			sum += xj * twiddle(j*k%n, n)
		}
		y[k] = sum
	}
	return
}

type transformVector struct {
	reference.Meta
	Params struct {
		X reference.Complex `json:"x"`
	} `json:"params"`
	Output struct {
		Y reference.Complex `json:"y"`
	} `json:"output"`
}

// The oracle must itself agree with SciPy before anything is tested against it.
// Measured: 0 to 2 ULP at small n, growing with n as its O(n) summations accumulate rounding,
// to 18 ULP, 4.0e-15, at n = 4096.
// Tolerance: about 10 times the worst.
func TestDFTMatchesSciPy(t *testing.T) {
	reference.Run(t, reference.SourceSciPy, "fft", func(t *testing.T, v transformVector) {
		x, err := v.Params.X.Values()
		require.NoError(t, err)
		want, err := v.Output.Y.Values()
		require.NoError(t, err)
		reference.AssertClose(t, "y", dft(x), want, reference.Tolerance{Rel: 1e-13, Scaled: 1e-13})
	})
}
