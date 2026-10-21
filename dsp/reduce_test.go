package dsp

import (
	"math"
	"slices"
	"testing"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

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
