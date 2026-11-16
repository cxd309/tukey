package dsp

import (
	"slices"
	"testing"

	"github.com/cxd309/tukey/internal/poly"
	"github.com/cxd309/tukey/internal/reference"
	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

// cascading the sections must give back the same filter as BA:
// multiplying each section's b (and a) polynomials together should reproduce b (and a)
func TestSOSCascadesToBA(t *testing.T) {
	for _, band := range []Band{Lowpass(0.3), Highpass(0.2), Bandpass(0.2, 0.5), Bandstop(0.2, 0.5)} {
		t.Run(band.String(), func(t *testing.T) {
			f, err := Butter(3, band) // odd order: exercises the first-order section
			require.NoError(t, err)
			b, a, err := f.BA()
			require.NoError(t, err)
			sos, err := f.SOS()
			require.NoError(t, err)

			cascadeB, cascadeA := []float64{1}, []float64{1}
			for _, s := range sos {
				cascadeB = poly.Mul(cascadeB, s[:3])
				cascadeA = poly.Mul(cascadeA, s[3:])
			}
			// padding sections carry a pole and zero at the origin, which appear as
			// trailing zeros; trim to BA's length before comparing
			tol := reference.Tolerance{Rel: 1e-13, Scaled: 1e-13}
			reference.AssertClose(t, "b", cascadeB[:len(b)], b, tol)
			reference.AssertClose(t, "a", cascadeA[:len(a)], a, tol)
		})
	}
}

func TestSOSRejectsUnpairedComplexRoots(t *testing.T) {
	_, err := ZPK{Poles: []complex128{complex(0.5, 0.5)}, Gain: 1}.SOS()
	assert.ErrorIs(t, err, ErrInvalidZPK)
}

func TestSOSDoesNotModifyZPK(t *testing.T) {
	f, err := Butter(4, Bandpass(0.2, 0.5))
	require.NoError(t, err)
	zeros, poles := slices.Clone(f.Zeros), slices.Clone(f.Poles)

	_, err = f.SOS()
	require.NoError(t, err)
	assert.Equal(t, zeros, f.Zeros)
	assert.Equal(t, poles, f.Poles)
}

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
