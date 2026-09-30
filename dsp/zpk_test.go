package dsp

import (
	"testing"

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
