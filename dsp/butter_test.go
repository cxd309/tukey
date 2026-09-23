package dsp

import (
	"math"
	"testing"

	"github.com/stretchr/testify/assert"
)

func TestButterRejectsInvalidInput(t *testing.T) {
	cases := []struct {
		name    string
		order   int
		band    Band
		wantErr string
	}{
		{"order 0", 0, Lowpass(0.3), "order must be >= 1"},
		{"cutoff 0", 2, Lowpass(0), "Lowpass(0): frequencies must be in (0, 1)"},
		{"cutoff 1", 2, Highpass(1), "Highpass(1): frequencies must be in (0, 1)"},
		{"cutoff NaN", 2, Lowpass(math.NaN()), "frequencies must be in (0, 1)"},
		{"high edge 1", 2, Bandstop(0.2, 1), "Bandstop(0.2, 1): frequencies must be in (0, 1)"},
		{"low >= high", 2, Bandpass(0.5, 0.2), "Bandpass(0.5, 0.2): low edge must be below high edge"},
		{"zero Band", 2, Band{}, "zero Band"},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			_, _, err := Butter(c.order, c.band)
			assert.ErrorContains(t, err, c.wantErr)
		})
	}
}
