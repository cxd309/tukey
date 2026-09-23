package reference

import (
	"math"
	"testing"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

func TestCompare(t *testing.T) {
	tol := Tolerance{Rel: 1e-12, Abs: 1e-12}
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

func TestCompareRejectsLengthMismatch(t *testing.T) {
	_, err := Compare([]float64{1}, []float64{1, 2}, Tolerance{})
	assert.ErrorContains(t, err, "length mismatch")
}
