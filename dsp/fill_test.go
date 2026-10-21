package dsp

import (
	"math"
	"slices"
	"testing"

	"github.com/stretchr/testify/assert"
)

var nan = math.NaN()

func TestFillLinear(t *testing.T) {
	cases := []struct {
		name string
		x    []float64
		want []float64
	}{
		{"nothing missing", []float64{1, 2, 4}, []float64{1, 2, 4}},
		{"gap", []float64{1, nan, 3}, []float64{1, 2, 3}},
		{"long gap", []float64{0, nan, nan, nan, 4}, []float64{0, 1, 2, 3, 4}},
		{"leading, extrapolated", []float64{nan, nan, 2, 3}, []float64{0, 1, 2, 3}},
		{"trailing, extrapolated", []float64{1, 2, nan, nan}, []float64{1, 2, 3, 4}},
		{"one known, unchanged", []float64{nan, 5, nan}, []float64{nan, 5, nan}},
		{"none known, unchanged", []float64{nan, nan}, []float64{nan, nan}},
		{"empty", []float64{}, []float64{}},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			assertEqualNaN(t, c.want, FillMissing(c.x, FillLinear()))
		})
	}
}

func TestFillConstant(t *testing.T) {
	assertEqualNaN(t, []float64{0, 1, 0, 2}, FillMissing([]float64{nan, 1, nan, 2}, FillConstant(0)))
	assertEqualNaN(t, []float64{7, 7}, FillMissing([]float64{nan, nan}, FillConstant(7)))
}

// The zero Fill is FillLinear.
func TestFillZeroIsLinear(t *testing.T) {
	x := []float64{nan, 1, nan, 3, nan}
	assertEqualNaN(t, FillMissing(x, FillLinear()), FillMissing(x, Fill{}))
}

// FillMissing returns a new slice: changing it, even when nothing was missing,
// doesn't change the input.
func TestFillMissingDoesNotModifyInput(t *testing.T) {
	for _, fill := range []Fill{FillLinear(), FillConstant(0)} {
		x := []float64{nan, 1, nan, 3}
		before := slices.Clone(x)
		y := FillMissing(x, fill)
		y[1] = 100
		assertEqualNaN(t, before, x)
	}
}

// assertEqualNaN asserts got equals want element by element, with NaN equal to
// NaN, which assert.Equal treats as different.
func assertEqualNaN(t *testing.T, want, got []float64) {
	t.Helper()
	if !assert.Len(t, got, len(want)) {
		return
	}
	for i := range want {
		if math.IsNaN(want[i]) {
			assert.True(t, math.IsNaN(got[i]), "got[%d] = %g, want NaN", i, got[i])
			continue
		}
		assert.Equal(t, want[i], got[i], "got[%d]", i)
	}
}
