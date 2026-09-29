package dsp

import (
	"testing"

	"github.com/cxd309/godsp/internal/reference"
	"github.com/stretchr/testify/assert"
)

// the cases use small integer roots, so results should be exact to rounding
var polyTolerance = reference.Tolerance{Rel: 1e-12, Scaled: 1e-12}

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
			reference.AssertClose(t, "coeffs", polyFromRoots(c.roots), c.want, polyTolerance)
		})
	}
}

func TestRealCoeffs(t *testing.T) {
	t.Run("clean values pass through", func(t *testing.T) {
		in := []complex128{complex(1, 0), complex(-2, 1e-15), complex(5, -1e-14)}
		want := []float64{1, -2, 5}
		reference.AssertClose(t, "coeffs", realCoeffs(in), want, polyTolerance)
	})

	t.Run("non-negligible imaginary part panics", func(t *testing.T) {
		assert.Panics(t, func() {
			realCoeffs([]complex128{complex(1, 0.5)})
		})
	})
}
