package dsp

import (
	"testing"

	"github.com/stretchr/testify/assert"
)

// Band.String formats a band as the constructor call that built it,
// which is what error messages show
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
