package reference

import (
	"math"
	"testing"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

func TestCompare(t *testing.T) {
	tol := Tolerance{Rel: 1e-12, Scaled: 1e-12}
	cases := []struct {
		name       string
		got, want  []float64
		mismatches int
	}{
		{"exact", []float64{1, 2}, []float64{1, 2}, 0},
		{"within tolerance", []float64{1 + 1e-15}, []float64{1}, 0},
		{"tiny values still checked", []float64{2e-14}, []float64{1e-14}, 1},
		{"NaN in got", []float64{math.NaN()}, []float64{1}, 1},
		{"NaN in both", []float64{math.NaN(), 1}, []float64{math.NaN(), 1}, 0},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			c, err := Compare(tc.got, tc.want, tol)
			require.NoError(t, err)
			assert.Len(t, c.Mismatches, tc.mismatches)
		})
	}
}

// the butter2_bp0.2-0.5_constant case: a bandpass rejecting a constant input of 3
// leaves only rounding noise, which relative to the ~0 output alone looks enormous
func TestCompareScaledFloor(t *testing.T) {
	tol := Tolerance{Rel: 1e-12, Scaled: 1e-12}
	got, want := []float64{5.8e-17}, []float64{3.85e-33}

	c, err := Compare(got, want, tol)
	require.NoError(t, err)
	assert.Len(t, c.Mismatches, 1, "relative to ~0 alone, rounding noise looks enormous")

	c, err = CompareScaled(got, want, tol, 3)
	require.NoError(t, err)
	assert.Empty(t, c.Mismatches, "relative to the input's scale, it's negligible")
}

func TestCompareRejectsLengthMismatch(t *testing.T) {
	_, err := Compare([]float64{1}, []float64{1, 2}, Tolerance{})
	assert.ErrorContains(t, err, "length mismatch")
}
