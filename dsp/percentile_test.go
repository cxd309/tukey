package dsp

import (
	"math"
	"slices"
	"testing"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

// The cases use x = 4, 1, 3, 2, unsorted, so 1, 2, 3, 4 sorted.
func TestPercentile(t *testing.T) {
	x := []float64{4, 1, 3, 2}
	cases := []struct {
		name   string
		method PercentileMethod
		p      float64
		want   float64
	}{
		// linear: position p/100·3
		{"linear minimum", PercentileLinear, 0, 1},
		{"linear maximum", PercentileLinear, 100, 4},
		{"linear median", PercentileLinear, 50, 2.5},
		{"linear quartile", PercentileLinear, 25, 1.75},
		// Hazen: position p/100·4 − 0.5, clamped to 0 to 3
		{"hazen median", PercentileHazen, 50, 2.5},
		{"hazen quartile", PercentileHazen, 25, 1.5},
		{"hazen clamped low", PercentileHazen, 10, 1},
		{"hazen clamped high", PercentileHazen, 95, 4},
		{"hazen minimum", PercentileHazen, 0, 1},
		{"hazen maximum", PercentileHazen, 100, 4},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			value, err := Percentile(x, c.p, PercentileBy(c.method))
			require.NoError(t, err)
			assert.Equal(t, c.want, value)
		})
	}
}

// The zero PercentileMethod, and so no option, is PercentileLinear.
func TestPercentileDefaultIsLinear(t *testing.T) {
	var method PercentileMethod
	assert.Equal(t, PercentileLinear, method)

	x := noise(100, 1)
	for _, p := range []float64{0, 12.5, 50, 99, 100} {
		byDefault, err := Percentile(x, p)
		require.NoError(t, err)
		linear, err := Percentile(x, p, PercentileBy(PercentileLinear))
		require.NoError(t, err)
		assert.Equal(t, linear, byDefault, "p %g", p)
	}
}

func TestPercentileOneSample(t *testing.T) {
	for _, method := range []PercentileMethod{PercentileLinear, PercentileHazen} {
		for _, p := range []float64{0, 30, 50, 100} {
			value, err := Percentile([]float64{7}, p, PercentileBy(method))
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

// As p rises the percentile never falls, and stays within x's range, for both
// methods: interpolating from the nearer neighbour keeps it monotonic.
func TestPercentileIsMonotonic(t *testing.T) {
	x := noise(101, 2)
	lowest, highest := slices.Min(x), slices.Max(x)
	for _, method := range []PercentileMethod{PercentileLinear, PercentileHazen} {
		previous := math.Inf(-1)
		for p := 0.0; p <= 100; p += 0.25 {
			value, err := Percentile(x, p, PercentileBy(method))
			require.NoError(t, err)
			assert.GreaterOrEqual(t, value, previous, "method %d, p %g", method, p)
			assert.GreaterOrEqual(t, value, lowest)
			assert.LessOrEqual(t, value, highest)
			previous = value
		}
	}
}

// Any NaN in x gives NaN, not an error, as a missing sample makes the
// percentile unknown. MATLAB's prctile ignores NaN instead.
func TestPercentileOfNaNIsNaN(t *testing.T) {
	for _, method := range []PercentileMethod{PercentileLinear, PercentileHazen} {
		value, err := Percentile([]float64{1, nan, 3, 2}, 50, PercentileBy(method))
		require.NoError(t, err)
		assert.True(t, math.IsNaN(value), "method %d", method)
	}
}

func TestPercentileRejectsInvalidInput(t *testing.T) {
	for _, x := range [][]float64{nil, {}} {
		_, err := Percentile(x, 50)
		assert.ErrorIs(t, err, ErrInvalidLength)
	}
	for _, p := range []float64{-0.001, 100.001, math.NaN(), math.Inf(1)} {
		_, err := Percentile([]float64{1, 2, 3}, p)
		assert.ErrorIs(t, err, ErrInvalidParameter, "p %g", p)
	}
	_, err := Percentile([]float64{1, 2, 3}, 50, PercentileBy(PercentileMethod(7)))
	assert.ErrorIs(t, err, ErrInvalidMode)
	assert.ErrorContains(t, err, "PercentileMethod 7")
}

// Percentile sorts a copy, leaving x in its order.
func TestPercentileDoesNotModifyInput(t *testing.T) {
	x := []float64{4, 1, 3, 2}
	before := slices.Clone(x)
	_, err := Percentile(x, 50)
	require.NoError(t, err)
	assert.Equal(t, before, x)
}
