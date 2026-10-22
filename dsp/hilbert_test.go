package dsp

import (
	"math"
	"testing"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

// The analytic signal's real part is x itself.
func TestHilbertRealPartIsX(t *testing.T) {
	for _, x := range [][]float64{
		{0.5},
		{1, -2},
		{3, 1, -4, 1, 5, -9, 2, 6},
		{3, 1, -4, 1, 5, -9, 2, 6, -5},
	} {
		analytic, err := Hilbert(x)
		require.NoError(t, err)
		require.Len(t, analytic, len(x))
		for i, a := range analytic {
			assert.InDelta(t, x[i], real(a), 1e-14, "real(analytic[%d]), n=%d", i, len(x))
		}
	}
}

// A cosine's Hilbert transform is the sine at the same frequency, so its
// analytic signal is a complex exponential, at odd and even lengths.
func TestHilbertCosineToSine(t *testing.T) {
	for _, n := range []int{15, 16} {
		x := make([]float64, n)
		for i := range x {
			x[i] = math.Cos(2 * math.Pi * 3 * float64(i) / float64(n))
		}
		analytic, err := Hilbert(x)
		require.NoError(t, err)
		for i, a := range analytic {
			want := math.Sin(2 * math.Pi * 3 * float64(i) / float64(n))
			assert.InDelta(t, want, imag(a), 1e-14, "imag(analytic[%d]), n=%d", i, n)
		}
	}
}

// DC and Nyquist are their own negative frequencies, so they're kept as they
// are, not doubled, and have no Hilbert transform.
func TestHilbertDCAndNyquist(t *testing.T) {
	cases := []struct {
		name string
		x    []float64
	}{
		{"DC", []float64{2, 2, 2, 2}},
		{"Nyquist", []float64{1, -1, 1, -1}},
		{"both, length 2", []float64{3, 1}},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			analytic, err := Hilbert(c.x)
			require.NoError(t, err)
			for i, a := range analytic {
				assert.InDelta(t, c.x[i], real(a), 1e-15, "real(analytic[%d])", i)
				assert.InDelta(t, 0, imag(a), 1e-15, "imag(analytic[%d])", i)
			}
		})
	}
}

func TestHilbertEmpty(t *testing.T) {
	_, err := Hilbert(nil)
	assert.ErrorIs(t, err, ErrInvalidLength)
}
