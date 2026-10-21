package window

import (
	"fmt"
	"testing"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

// windows is every window, with any parameters fixed
// for the tests every window should pass
var windows = map[string]func(n int, opts ...Option) ([]float64, error){
	"Hann":           Hann,
	"Hamming":        Hamming,
	"Blackman":       Blackman,
	"BlackmanHarris": BlackmanHarris,
	"Nuttall":        Nuttall,
	"FlatTop":        FlatTop,
	"GeneralHamming(0.75)": func(n int, opts ...Option) ([]float64, error) {
		return GeneralHamming(n, 0.75, opts...)
	},
	"GeneralCosine(HFT90D)": func(n int, opts ...Option) ([]float64, error) {
		return GeneralCosine(n, []float64{1, 1.942604, 1.340318, 0.440811, 0.043097}, opts...)
	},
	"Boxcar":   Boxcar,
	"Triang":   Triang,
	"Bartlett": Bartlett,
	"Parzen":   Parzen,
	"Bohman":   Bohman,
	"Barthann": Barthann,
	"Cosine":   Cosine,
	"Lanczos":  Lanczos,
	"Hanning":  Hanning,
	"Tukey(0.25)": func(n int, opts ...Option) ([]float64, error) {
		return Tukey(n, 0.25, opts...)
	},
	"Tukey(0.5)": func(n int, opts ...Option) ([]float64, error) {
		return Tukey(n, 0.5, opts...)
	},
	"Gaussian(2.5)": func(n int, opts ...Option) ([]float64, error) {
		return Gaussian(n, 2.5, opts...)
	},
	"GeneralGaussian(1.5, 7)": func(n int, opts ...Option) ([]float64, error) {
		return GeneralGaussian(n, 1.5, 7, opts...)
	},
	"Exponential(3)": func(n int, opts ...Option) ([]float64, error) {
		return Exponential(n, 3, opts...)
	},
	"Kaiser(8.6)": func(n int, opts ...Option) ([]float64, error) {
		return Kaiser(n, 8.6, opts...)
	},
	"Taylor(4, 30)": func(n int, opts ...Option) ([]float64, error) {
		return Taylor(n, 4, 30, opts...)
	},
	"Taylor(9, 50, Unnormalised)": func(n int, opts ...Option) ([]float64, error) {
		return Taylor(n, 9, 50, append(opts, Unnormalised())...)
	},
	"Chebwin(50)": func(n int, opts ...Option) ([]float64, error) {
		return Chebwin(n, 50, opts...)
	},
}

func TestRejectsNegativeLength(t *testing.T) {
	for name, window := range windows {
		t.Run(name, func(t *testing.T) {
			_, err := window(-1)
			assert.ErrorIs(t, err, ErrInvalidLength)
		})
	}
}

// n = 0 and n = 1 are fixed for every window, before any formula, as SciPy
func TestShortLengths(t *testing.T) {
	for name, window := range windows {
		for _, opts := range [][]Option{nil, {Periodic()}} {
			t.Run(fmt.Sprintf("%s/%d options", name, len(opts)), func(t *testing.T) {
				w, err := window(0, opts...)
				require.NoError(t, err)
				assert.Empty(t, w)

				w, err = window(1, opts...)
				require.NoError(t, err)
				assert.Equal(t, []float64{1}, w)
			})
		}
	}
}

// symmetric windows are mirrored, so symmetric exactly, not just to within rounding
func TestSymmetricIsExact(t *testing.T) {
	for name, window := range windows {
		t.Run(name, func(t *testing.T) {
			for n := 2; n <= 20; n++ {
				w, err := window(n)
				require.NoError(t, err)
				for k := range w {
					assert.Equal(t, w[k], w[n-1-k], "n = %d, k = %d", n, k)
				}
			}
		})
	}
}

// periodic is the symmetric window one longer, without its last point
func TestPeriodicIsSymmetricWithoutLastPoint(t *testing.T) {
	for name, window := range windows {
		if name == "Hanning" {
			continue // MATLAB defines periodic Hanning as Hann's periodic window: see TestHanning
		}
		t.Run(name, func(t *testing.T) {
			for n := 2; n <= 20; n++ {
				periodic, err := window(n, Periodic())
				require.NoError(t, err)
				symmetric, err := window(n + 1)
				require.NoError(t, err)
				assert.Equal(t, symmetric[:n], periodic, "n = %d", n)
			}
		})
	}
}

// assertInvalidParameter checks that each call returns ErrInvalidParameter
func assertInvalidParameter(t *testing.T, calls map[string]func() ([]float64, error)) {
	t.Helper()
	for name, call := range calls {
		t.Run(name, func(t *testing.T) {
			_, err := call()
			assert.ErrorIs(t, err, ErrInvalidParameter)
		})
	}
}
