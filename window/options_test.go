package window

import (
	"testing"

	"github.com/stretchr/testify/assert"
)

// each window accepts only the options that apply to it
func TestRejectsOptionThatDoesNotApply(t *testing.T) {
	cases := []struct {
		name    string
		call    func() ([]float64, error)
		message string
	}{
		{"Center", func() ([]float64, error) { return Hann(8, Center(2)) }, "Center doesn't apply to Hann"},
		{"Unnormalised", func() ([]float64, error) { return Hann(8, Unnormalised()) }, "Unnormalised doesn't apply to Hann"},
		{"Periodic", func() ([]float64, error) { return KaiserBesselDerived(8, 4, Periodic()) }, "Periodic doesn't apply to KaiserBesselDerived"},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			_, err := c.call()
			assert.ErrorIs(t, err, ErrInvalidParameter)
			assert.ErrorContains(t, err, c.message)
		})
	}
}
