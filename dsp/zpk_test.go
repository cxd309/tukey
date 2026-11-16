package dsp

import (
	"slices"
	"testing"

	"github.com/cxd309/tukey/internal/poly"
	"github.com/cxd309/tukey/internal/reference"
	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

// a caller-built ZPK with a conjugate pair expands to real coefficients:
// zeros +-j and poles 0.5+-0.5j give b= [1, 0, 1], a= [1, -1, 0.5]
func TestZPKBAConjugatePairs(t *testing.T) {
	f := ZPK{
		Zeros: []complex128{1i, -1i},
		Poles: []complex128{complex(0.5, 0.5), complex(0.5, -0.5)},
		Gain:  1,
	}
	b, a, err := f.BA()
	require.NoError(t, err)
	assert.Equal(t, []float64{1, 0, 1}, b)
	assert.Equal(t, []float64{1, -1, 0.5}, a)
}

func TestZPKBARejectsUnpairedComplexRoots(t *testing.T) {
	cases := []struct {
		name string
		f    ZPK
	}{
		{"unpaired zero", ZPK{Zeros: []complex128{1i}, Gain: 1}},
		{"unpaired pole", ZPK{Poles: []complex128{complex(0.5, 0.5)}, Gain: 1}},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			_, _, err := c.f.BA()
			assert.ErrorIs(t, err, ErrInvalidZPK)
		})
	}
}

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
