package fft

import (
	"testing"

	"github.com/stretchr/testify/assert"
)

func TestFreq(t *testing.T) {
	// even n: Nyquist (0.5) is listed as negative
	assert.Equal(t, []float64{0, 0.125, 0.25, 0.375, -0.5, -0.375, -0.25, -0.125}, Freq(8, 1))
	// odd n: no Nyquist bin, so positives and negatives balance
	assert.Equal(t, []float64{0, 0.2, 0.4, -0.4, -0.2}, Freq(5, 1))
	// d = 1/fs gives Hz: 4 samples at 1 kHz
	assert.Equal(t, []float64{0, 250, -500, -250}, Freq(4, 0.001))
	assert.Equal(t, []float64{0}, Freq(1, 1))
	assert.Empty(t, Freq(0, 1))
}

func TestRealFreq(t *testing.T) {
	// even n: ends at +Nyquist
	assert.Equal(t, []float64{0, 0.125, 0.25, 0.375, 0.5}, RealFreq(8, 1))
	// odd n: stops short of Nyquist
	assert.Equal(t, []float64{0, 0.2, 0.4}, RealFreq(5, 1))
	assert.Equal(t, []float64{0, 250, 500}, RealFreq(4, 0.001))
	assert.Equal(t, []float64{0}, RealFreq(1, 1))
	assert.Empty(t, RealFreq(0, 1))
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

// Shift and InverseShift return new slices, never the caller's
func TestShiftDoesNotAlias(t *testing.T) {
	x := []int{0, 1, 2, 3}
	Shift(x)[0] = 99
	InverseShift(x)[0] = 99
	assert.Equal(t, []int{0, 1, 2, 3}, x)
}
