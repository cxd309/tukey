package window

import (
	"fmt"
	"math"
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
			continue // MATLAB defines periodic Hanning as Hann's periodic window
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

// the general windows reduce to the named ones exactly: they share one code path
func TestSpecialCasesAreExact(t *testing.T) {
	for n := 0; n <= 20; n++ {
		hann, err := Hann(n)
		require.NoError(t, err)
		hamming, err := Hamming(n)
		require.NoError(t, err)

		w, err := GeneralHamming(n, 0.5)
		require.NoError(t, err)
		assert.Equal(t, hann, w, "GeneralHamming(%d, 0.5)", n)

		w, err = GeneralHamming(n, 0.54)
		require.NoError(t, err)
		assert.Equal(t, hamming, w, "GeneralHamming(%d, 0.54)", n)

		w, err = GeneralCosine(n, []float64{0.5, 0.5})
		require.NoError(t, err)
		assert.Equal(t, hann, w, "GeneralCosine(%d, [0.5, 0.5])", n)
	}
}

func TestHannSmall(t *testing.T) {
	w, err := Hann(3)
	require.NoError(t, err)
	assert.Equal(t, []float64{0, 1, 0}, w)

	w, err = Hann(4, Periodic())
	require.NoError(t, err)
	assert.InDeltaSlice(t, []float64{0, 0.5, 1, 0.5}, w, 1e-16)
}

func TestGeneralCosineNeedsCoefficients(t *testing.T) {
	_, err := GeneralCosine(8, nil)
	assert.ErrorIs(t, err, ErrInvalidParameter)
}

func TestRejectsOptionThatDoesNotApply(t *testing.T) {
	_, err := newConfig("Test", []Option{Periodic()}, optionSet{})
	assert.ErrorIs(t, err, ErrInvalidParameter)
	assert.ErrorContains(t, err, "Periodic doesn't apply to Test")
}

// Tukey's alpha has no domain: at or below 0 it's Boxcar, at or above 1 it's Hann
func TestTukeyLimitsAreExact(t *testing.T) {
	for n := 0; n <= 20; n++ {
		boxcar, err := Boxcar(n)
		require.NoError(t, err)
		hann, err := Hann(n)
		require.NoError(t, err)

		for _, alpha := range []float64{-1, 0} {
			w, err := Tukey(n, alpha)
			require.NoError(t, err)
			assert.Equal(t, boxcar, w, "Tukey(%d, %g)", n, alpha)
		}
		for _, alpha := range []float64{1, 2} {
			w, err := Tukey(n, alpha)
			require.NoError(t, err)
			assert.Equal(t, hann, w, "Tukey(%d, %g)", n, alpha)
		}
	}
}

func TestTukeyRejectsNaN(t *testing.T) {
	_, err := Tukey(8, math.NaN())
	assert.ErrorIs(t, err, ErrInvalidParameter)
}

// Hanning is Hann without its zero endpoints; periodic, it's Hann's periodic window
func TestHanning(t *testing.T) {
	for n := 0; n <= 20; n++ {
		hann, err := Hann(n + 2)
		require.NoError(t, err)
		w, err := Hanning(n)
		require.NoError(t, err)
		assert.Equal(t, hann[1:n+1], w, "Hanning(%d)", n)

		hann, err = Hann(n, Periodic())
		require.NoError(t, err)
		w, err = Hanning(n, Periodic())
		require.NoError(t, err)
		assert.Equal(t, hann, w, "Hanning(%d, Periodic())", n)
	}
}

func TestTriangularSmall(t *testing.T) {
	w, err := Triang(3)
	require.NoError(t, err)
	assert.Equal(t, []float64{0.5, 1, 0.5}, w)

	w, err = Triang(4)
	require.NoError(t, err)
	assert.Equal(t, []float64{0.25, 0.75, 0.75, 0.25}, w)

	w, err = Bartlett(3)
	require.NoError(t, err)
	assert.Equal(t, []float64{0, 1, 0}, w)

	w, err = Bohman(5)
	require.NoError(t, err)
	assert.Equal(t, 0.0, w[0])
	assert.Equal(t, 0.0, w[4])
}
