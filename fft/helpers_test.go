package fft

import (
	"math"
	"testing"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

func TestFreq(t *testing.T) {
	cases := []struct {
		name string
		n    int
		d    float64
		want []float64
	}{
		{"even n lists Nyquist as negative", 8, 1, []float64{0, 0.125, 0.25, 0.375, -0.5, -0.375, -0.25, -0.125}},
		{"odd n has no Nyquist bin", 5, 1, []float64{0, 0.2, 0.4, -0.4, -0.2}},
		{"d = 1/fs gives Hz, 4 samples at 1 kHz", 4, 0.001, []float64{0, 250, -500, -250}},
		{"one bin", 1, 1, []float64{0}},
		{"no bins", 0, 1, []float64{}},
	}
	for _, c := range cases {
		f, err := Freq(c.n, c.d)
		require.NoError(t, err, c.name)
		assert.Equal(t, c.want, f, c.name)
	}
}

func TestRealFreq(t *testing.T) {
	cases := []struct {
		name string
		n    int
		d    float64
		want []float64
	}{
		{"even n ends at Nyquist", 8, 1, []float64{0, 0.125, 0.25, 0.375, 0.5}},
		{"odd n stops short of Nyquist", 5, 1, []float64{0, 0.2, 0.4}},
		{"d = 1/fs gives Hz", 4, 0.001, []float64{0, 250, 500}},
		{"one bin", 1, 1, []float64{0}},
		{"no bins", 0, 1, []float64{}},
	}
	for _, c := range cases {
		f, err := RealFreq(c.n, c.d)
		require.NoError(t, err, c.name)
		assert.Equal(t, c.want, f, c.name)
	}
}

func TestFreqRejectsInvalidInput(t *testing.T) {
	for _, freq := range []func(int, float64) ([]float64, error){Freq, RealFreq} {
		_, err := freq(-1, 1)
		assert.ErrorIs(t, err, ErrInvalidLength)
		for _, d := range []float64{0, -1, math.NaN(), math.Inf(1)} {
			_, err = freq(4, d)
			assert.ErrorIs(t, err, ErrInvalidParameter, "d %g", d)
		}
	}
}

func TestShift(t *testing.T) {
	assert.Equal(t, []int{2, 3, 0, 1}, Shift([]int{0, 1, 2, 3}))
	assert.Equal(t, []int{3, 4, 0, 1, 2}, Shift([]int{0, 1, 2, 3, 4}))
	assert.Equal(t, []int{0}, Shift([]int{0}))
	assert.Empty(t, Shift([]int{}))
}

func TestInverseShift(t *testing.T) {
	assert.Equal(t, []int{2, 3, 0, 1}, InverseShift([]int{0, 1, 2, 3}))       // even: same as Shift
	assert.Equal(t, []int{2, 3, 4, 0, 1}, InverseShift([]int{0, 1, 2, 3, 4})) // odd: differs
	assert.Empty(t, InverseShift([]int{}))
}

// An empty or nil x gives an empty, non-nil result.
func TestShiftEmptyIsNotNil(t *testing.T) {
	for _, x := range [][]int{nil, {}} {
		assert.NotNil(t, Shift(x))
		assert.NotNil(t, InverseShift(x))
	}
}

// Shift and InverseShift return new slices, never the caller's
func TestShiftDoesNotAlias(t *testing.T) {
	x := []int{0, 1, 2, 3}
	Shift(x)[0] = 99
	InverseShift(x)[0] = 99
	assert.Equal(t, []int{0, 1, 2, 3}, x)
}
