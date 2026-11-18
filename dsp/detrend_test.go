package dsp

import (
	"slices"
	"testing"

	"github.com/stretchr/testify/assert"
)

// A straight line with exactly representable samples detrends to exactly zero.
// Every intermediate is exact:
// the mean is 4, the centred times are whole numbers, and the fitted slope is exactly 0.5.
func TestDetrendLinearRemovesALineExactly(t *testing.T) {
	x := make([]float64, 9)
	for i := range x {
		x[i] = 2 + 0.5*float64(i)
	}
	y := Detrend(x, DetrendLinear())
	assert.Equal(t, make([]float64, 9), y)
}

func TestDetrendConstantRemovesTheMean(t *testing.T) {
	y := Detrend([]float64{1, 2, 3, 6}, DetrendConstant())
	assert.Equal(t, []float64{-2, -1, 0, 3}, y)
}

// The zero DetrendMode is DetrendLinear(), same as SciPy's and MATLAB's default.
func TestDetrendZeroModeIsLinear(t *testing.T) {
	var mode DetrendMode
	assert.Equal(t, DetrendLinear(), mode)
}

// detrendModes is both modes, named for subtests.
var detrendModes = []struct {
	name string
	mode DetrendMode
}{
	{"linear", DetrendLinear()},
	{"constant", DetrendConstant()},
}

func TestDetrendShortSignals(t *testing.T) {
	for _, m := range detrendModes {
		t.Run(m.name, func(t *testing.T) {
			y := Detrend(nil, m.mode)
			assert.Empty(t, y)

			y = Detrend([]float64{5}, m.mode)
			assert.Equal(t, []float64{0}, y, "a single sample is its own mean, with no slope")
		})
	}
}

func TestDetrendDoesNotModifyInput(t *testing.T) {
	x := []float64{1, 4, 2, 8}
	before := slices.Clone(x)
	for _, m := range detrendModes {
		Detrend(x, m.mode)
		assert.Equal(t, before, x, m.name)
	}
}
