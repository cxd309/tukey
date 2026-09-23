package dsp

import (
	"math"
	"testing"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

func TestPolyFromRoots(t *testing.T) {
	cases := []struct {
		name  string
		roots []complex128
		want  []complex128 // highest degree first
	}{
		{"no roots", nil, []complex128{1}},
		{"single real root", []complex128{1}, []complex128{1, -1}},
		{"two real roots (x-2)(x-3) = x^2 -5x +6", []complex128{2, 3}, []complex128{1, -5, 6}},
		{
			"conjugate pair (x-(1+2i))(x-(1-2i)) = x^2 -2x +5",
			[]complex128{complex(1, 2), complex(1, -2)},
			[]complex128{1, -2, 5},
		},
	}

	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			got := polyFromRoots(c.roots)
			require.Len(t, got, len(c.want))
			for i := range c.want {
				assert.LessOrEqual(t, math.Hypot(real(got[i]-c.want[i]), imag(got[i]-c.want[i])), 1e-9,
					"coeff %d: got %v, want %v", i, got[i], c.want[i])
			}
		})
	}
}

func TestRealCoeffs(t *testing.T) {
	t.Run("clean values pass through", func(t *testing.T) {
		in := []complex128{complex(1, 0), complex(-2, 1e-15), complex(5, -1e-14)}
		want := []float64{1, -2, 5}
		got := realCoeffs(in)
		assert.InDeltaSlice(t, want, got, 1e-9)
	})

	t.Run("non-negligible imaginary part panics", func(t *testing.T) {
		assert.Panics(t, func() {
			realCoeffs([]complex128{complex(1, 0.5)})
		})
	})
}
