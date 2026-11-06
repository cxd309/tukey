package dsp

import (
	"math"
	"slices"
	"testing"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

// mixedMagnitudes is 1e8, then 20,000 samples of 0.001, then −1e8: its exact
// total is 20, but a plain sum loses the small values' low bits against 1e8.
func mixedMagnitudes() (x []float64) {
	x = []float64{1e8}
	for range 20_000 {
		x = append(x, 1e-3)
	}
	x = append(x, -1e8)
	return
}

func TestSum(t *testing.T) {
	cases := []struct {
		name string
		x    []float64
		want float64
	}{
		{"ramp", ramp(100), 5050},
		{"one sample", []float64{-7}, -7},
		{"empty", []float64{}, 0},
		{"nil", nil, 0},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			assert.Equal(t, c.want, Sum(c.x))
		})
	}
}

// Sum is compensated, so small values among large ones keep their low bits,
// where a plain sum drifts.
func TestSumKeepsSmallValuesAmongLarge(t *testing.T) {
	x := mixedMagnitudes()
	var plain float64
	for _, v := range x {
		plain += v
	}
	assert.NotEqual(t, 20.0, plain, "a plain sum should drift here, or the case tests nothing")
	assert.Equal(t, 20.0, Sum(x))
}

// 20,000 tenths sum to the double nearest 2000, which is 2000, where a plain
// sum's error grows with each addition.
func TestSumOfManyTenths(t *testing.T) {
	x := make([]float64, 20_000)
	for i := range x {
		x[i] = 0.1
	}
	assert.Equal(t, 2000.0, Sum(x))
}

func TestSumOfNaNIsNaN(t *testing.T) {
	assert.True(t, math.IsNaN(Sum([]float64{1, nan, 3})))
}

// An infinite sample or an overflowing total gives an infinite Sum, as a plain
// sum does, and opposite infinities give NaN.
func TestSumOfInfinity(t *testing.T) {
	assert.Equal(t, math.Inf(1), Sum([]float64{1, math.Inf(1), 2}))
	assert.Equal(t, math.Inf(-1), Sum([]float64{math.Inf(-1), 1}))
	assert.Equal(t, math.Inf(1), Sum([]float64{math.MaxFloat64, math.MaxFloat64}))
	assert.True(t, math.IsNaN(Sum([]float64{math.Inf(1), math.Inf(-1)})))
}

func TestMean(t *testing.T) {
	cases := []struct {
		name string
		x    []float64
		want float64
	}{
		{"whole", []float64{1, 2, 3, 6}, 3},
		{"fraction", []float64{1, 2}, 1.5},
		{"one sample", []float64{-7}, -7},
		{"mixed magnitudes", mixedMagnitudes(), 20.0 / 20_002},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			m, err := Mean(c.x)
			require.NoError(t, err)
			assert.Equal(t, c.want, m)
		})
	}
}

func TestMeanRejectsEmpty(t *testing.T) {
	for _, x := range [][]float64{nil, {}} {
		_, err := Mean(x)
		assert.ErrorIs(t, err, ErrInvalidLength)
	}
}

func TestRMS(t *testing.T) {
	cases := []struct {
		name string
		x    []float64
		want float64
	}{
		{"constant", []float64{-2, -2, -2}, 2},
		{"3-4", []float64{3, 4}, math.Sqrt(12.5)},
		{"one sample", []float64{-7}, 7},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			rms, err := RMS(c.x)
			require.NoError(t, err)
			assert.Equal(t, c.want, rms)
		})
	}
}

func TestRMSOfNaNIsNaN(t *testing.T) {
	rms, err := RMS([]float64{1, math.NaN()})
	require.NoError(t, err)
	assert.True(t, math.IsNaN(rms))
}

func TestRMSRejectsEmpty(t *testing.T) {
	_, err := RMS(nil)
	assert.ErrorIs(t, err, ErrInvalidLength)
}

func TestRMSDoesNotModifyInput(t *testing.T) {
	x := []float64{1, -2, 3}
	before := slices.Clone(x)
	_, err := RMS(x)
	require.NoError(t, err)
	assert.Equal(t, before, x)
}
