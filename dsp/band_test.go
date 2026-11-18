package dsp

import (
	"testing"

	"github.com/stretchr/testify/assert"
)

// Band.String formats a band as the constructor call that built it, which is what error messages show.
func TestBandString(t *testing.T) {
	cases := []struct {
		band Band
		want string
	}{
		{Lowpass(0.3), "Lowpass(0.3)"},
		{Highpass(0.7), "Highpass(0.7)"},
		{Bandpass(0.2, 0.5), "Bandpass(0.2, 0.5)"},
		{Bandstop(0.2, 0.5), "Bandstop(0.2, 0.5)"},
		{Band{}, "Band{}"},
	}
	for _, c := range cases {
		t.Run(c.want, func(t *testing.T) {
			assert.Equal(t, c.want, c.band.String())
		})
	}
}

// Edges is a copy, so changing it leaves the band as it was, and the zero Band's is empty rather than nil.
func TestBandEdges(t *testing.T) {
	band := Bandpass(0.2, 0.5)
	edges := band.Edges()
	assert.Equal(t, []float64{0.2, 0.5}, edges)
	edges[0] = 0.9
	assert.Equal(t, []float64{0.2, 0.5}, band.Edges())
	assert.Equal(t, []float64{0.3}, Lowpass(0.3).Edges())

	none := Band{}.Edges()
	assert.NotNil(t, none)
	assert.Empty(t, none)
}
