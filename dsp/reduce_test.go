package dsp

import (
	"math"
	"slices"
	"testing"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

// mixedMagnitudes is 1e8, then 20,000 samples of 0.001, then -1e8:
// its exact total is 20, but a plain sum loses the small values' low bits against 1e8.
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

// Sum is compensated, so small values among large ones keep their low bits, where a plain sum drifts.
func TestSumKeepsSmallValuesAmongLarge(t *testing.T) {
	x := mixedMagnitudes()
	var plain float64
	for _, v := range x {
		plain += v
	}
	assert.NotEqual(t, 20.0, plain, "a plain sum should drift here, or the case tests nothing")
	assert.Equal(t, 20.0, Sum(x))
}

// 20,000 tenths sum to the double nearest 2000, which is 2000,
// where a plain sum's error grows with each addition.
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

// An infinite sample or an overflowing total gives an infinite Sum, as a plain sum does,
// and opposite infinities give NaN.
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
		assert.ErrorIs(t, err, ErrSignalTooShort)
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
	assert.ErrorIs(t, err, ErrSignalTooShort)
}

func TestRMSDoesNotModifyInput(t *testing.T) {
	x := []float64{1, -2, 3}
	before := slices.Clone(x)
	_, err := RMS(x)
	require.NoError(t, err)
	assert.Equal(t, before, x)
}

// The cases use x = 4, 1, 3, 2, unsorted, so 1, 2, 3, 4 sorted.
func TestPercentile(t *testing.T) {
	x := []float64{4, 1, 3, 2}
	cases := []struct {
		name   string
		method PercentileOption
		p      float64
		want   float64
	}{
		// linear: position p/100·3
		{"linear minimum", PercentileLinear(), 0, 1},
		{"linear maximum", PercentileLinear(), 100, 4},
		{"linear median", PercentileLinear(), 50, 2.5},
		{"linear quartile", PercentileLinear(), 25, 1.75},
		// Hazen: position p/100·4 - 0.5, clamped to 0 to 3
		{"hazen median", PercentileHazen(), 50, 2.5},
		{"hazen quartile", PercentileHazen(), 25, 1.5},
		{"hazen clamped low", PercentileHazen(), 10, 1},
		{"hazen clamped high", PercentileHazen(), 95, 4},
		{"hazen minimum", PercentileHazen(), 0, 1},
		{"hazen maximum", PercentileHazen(), 100, 4},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			value, err := Percentile(x, c.p, c.method)
			require.NoError(t, err)
			assert.Equal(t, c.want, value)
		})
	}
}

// With no option, Percentile is PercentileLinear.
func TestPercentileDefaultIsLinear(t *testing.T) {
	x := noise(100, 1)
	for _, p := range []float64{0, 12.5, 50, 99, 100} {
		byDefault, err := Percentile(x, p)
		require.NoError(t, err)
		linear, err := Percentile(x, p, PercentileLinear())
		require.NoError(t, err)
		assert.Equal(t, linear, byDefault, "p %g", p)
	}
}

func TestPercentileOneSample(t *testing.T) {
	for _, method := range []PercentileOption{PercentileLinear(), PercentileHazen()} {
		for _, p := range []float64{0, 30, 50, 100} {
			value, err := Percentile([]float64{7}, p, method)
			require.NoError(t, err)
			assert.Equal(t, 7.0, value, "method %d, p %g", method, p)
		}
	}
}

// Between equal samples the percentile is that value exactly.
func TestPercentileTies(t *testing.T) {
	value, err := Percentile([]float64{2, 5, 2, 2}, 50)
	require.NoError(t, err)
	assert.Equal(t, 2.0, value)
}

// As p rises the percentile never falls, and stays within x's range, for both methods:
// interpolating from the nearer neighbour keeps it monotonic.
func TestPercentileIsMonotonic(t *testing.T) {
	x := noise(101, 2)
	lowest, highest := slices.Min(x), slices.Max(x)
	for _, method := range []PercentileOption{PercentileLinear(), PercentileHazen()} {
		previous := math.Inf(-1)
		for p := 0.0; p <= 100; p += 0.25 {
			value, err := Percentile(x, p, method)
			require.NoError(t, err)
			assert.GreaterOrEqual(t, value, previous, "method %d, p %g", method, p)
			assert.GreaterOrEqual(t, value, lowest)
			assert.LessOrEqual(t, value, highest)
			previous = value
		}
	}
}

// Any NaN in x gives NaN, not an error, as a missing sample makes the percentile unknown.
// MATLAB's prctile ignores NaN instead.
func TestPercentileOfNaNIsNaN(t *testing.T) {
	for _, method := range []PercentileOption{PercentileLinear(), PercentileHazen()} {
		value, err := Percentile([]float64{1, nan, 3, 2}, 50, method)
		require.NoError(t, err)
		assert.True(t, math.IsNaN(value), "method %d", method)
	}
}

func TestPercentileRejectsInvalidInput(t *testing.T) {
	for _, x := range [][]float64{nil, {}} {
		_, err := Percentile(x, 50)
		assert.ErrorIs(t, err, ErrSignalTooShort)
	}
	for _, p := range []float64{-0.001, 100.001, math.NaN(), math.Inf(1)} {
		_, err := Percentile([]float64{1, 2, 3}, p)
		assert.ErrorIs(t, err, ErrInvalidParameter, "p %g", p)
	}
}

// Percentile sorts a copy, leaving x in its order.
func TestPercentileDoesNotModifyInput(t *testing.T) {
	x := []float64{4, 1, 3, 2}
	before := slices.Clone(x)
	_, err := Percentile(x, 50)
	require.NoError(t, err)
	assert.Equal(t, before, x)
}

// Next to an infinite sample the percentile is that infinity,
// where NumPy's arithmetic gives NaN from Inf - Inf or 0·Inf.
func TestPercentileOfInfinity(t *testing.T) {
	inf := math.Inf(1)
	cases := []struct {
		name string
		x    []float64
		p    float64
		want float64
	}{
		{"minimum of -Inf", []float64{-inf, 1, 2, 3}, 0, -inf},
		{"maximum of +Inf", []float64{1, 2, inf}, 100, inf},
		{"between 1 and +Inf", []float64{1, inf}, 50, inf},
		{"between -Inf and 1", []float64{-inf, 1}, 50, -inf},
		{"a finite sample beside infinities", []float64{-inf, 2, inf}, 50, 2},
	}
	for _, method := range []PercentileOption{PercentileLinear(), PercentileHazen()} {
		for _, c := range cases {
			value, err := Percentile(c.x, c.p, method)
			require.NoError(t, err)
			assert.Equal(t, c.want, value, c.name)
		}
	}
	// halfway between opposite infinities has no value
	value, err := Percentile([]float64{-inf, inf}, 50)
	require.NoError(t, err)
	assert.True(t, math.IsNaN(value))
}
