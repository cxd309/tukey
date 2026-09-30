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
		wantErr error
	}{
		{"order 0", 0, Lowpass(0.3), ErrInvalidOrder},
		{"cutoff 0", 2, Lowpass(0), ErrInvalidBand},
		{"cutoff 1", 2, Highpass(1), ErrInvalidBand},
		{"cutoff NaN", 2, Lowpass(math.NaN()), ErrInvalidBand},
		{"high edge 1", 2, Bandstop(0.2, 1), ErrInvalidBand},
		{"low >= high", 2, Bandpass(0.5, 0.2), ErrInvalidBand},
		{"zero Band", 2, Band{}, ErrInvalidBand},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			_, err := Butter(c.order, c.band)
			assert.ErrorIs(t, err, c.wantErr)
		})
	}
}
