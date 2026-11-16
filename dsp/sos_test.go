package dsp

import (
	"testing"

	"github.com/stretchr/testify/assert"
)

func TestSOSNtaps(t *testing.T) {
	cases := []struct {
		name  string
		order int
		want  int
	}{
		{"even order: 2 sections", 4, 5},
		{"odd order: 2 sections, one first-order", 3, 4},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			assert.Equal(t, c.want, testSOS(t, c.order, Lowpass(0.3)).ntaps())
		})
	}
}

// SciPy's tap count and MATLAB's order differ for an odd-order Butterworth,
// whose real pole shares a section with two zeros: first-order to SciPy, order
// 2 to MATLAB. They agree for an even order.
func TestSOSSize(t *testing.T) {
	cases := []struct {
		name string
		sos  SOS
		want filterSize
	}{
		{"butter 3", testSOS(t, 3, Lowpass(0.3)), filterSize{ntaps: 4, order: 4}},
		{"butter 4", testSOS(t, 4, Lowpass(0.3)), filterSize{ntaps: 5, order: 4}},
		{"butter 2 bandstop", testSOS(t, 2, Bandstop(0.2, 0.5)), filterSize{ntaps: 5, order: 4}},
		// b2 = 0 in one section and a2 = 0 in the other: SciPy counts one
		// first-order section, MATLAB none
		{"mixed zeros", SOS{{0.2, 0.3, 0, 1, -0.5, 0.2}, {0.4, 0.1, 0.05, 1, -0.3, 0}}, filterSize{ntaps: 4, order: 4}},
		{"gain section", SOS{{2, 0, 0, 1, 0, 0}, {1, 1, 0, 1, -0.8, 0.4}}, filterSize{ntaps: 4, order: 2}},
		{"gain only", SOS{{2, 0, 0, 1, 0, 0}}, filterSize{ntaps: 2, order: 0}},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			assert.Equal(t, c.want, c.sos.size())
		})
	}
}
