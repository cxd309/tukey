package poly

import (
	"cmp"
	"math/cmplx"
	"slices"
	"testing"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

// drain takes every section from pr, in order.
func drain(pr Pairing, sections int) (zeros, poles [][]complex128) {
	for range sections {
		z, p := pr.Next()
		zeros, poles = append(zeros, z), append(poles, p)
	}
	return
}

// sortedRoots is roots sorted by real part then imaginary part, to compare as sets.
func sortedRoots(roots []complex128) (sorted []complex128) {
	sorted = slices.Clone(roots)
	slices.SortFunc(sorted, func(a, b complex128) int {
		if c := cmp.Compare(real(a), real(b)); c != 0 {
			return c
		}
		return cmp.Compare(imag(a), imag(b))
	})
	return
}

// Unequal counts are balanced with roots at the origin, and an odd count gets one more of each,
// so every section has two zeros and two poles, and the sections hold every root.
func TestPairingBalancesCounts(t *testing.T) {
	zeros := []complex128{-1}
	poles := []complex128{complex(0.5, 0.5), complex(0.5, -0.5), 0.9}
	pr, sections, err := NewPairing(zeros, poles)
	require.NoError(t, err)
	assert.Equal(t, 2, sections)
	gotZeros, gotPoles := drain(pr, sections)
	var allZeros, allPoles []complex128
	for i := range sections {
		assert.Len(t, gotZeros[i], 2, "section %d's zeros", i)
		assert.Len(t, gotPoles[i], 2, "section %d's poles", i)
		allZeros, allPoles = append(allZeros, gotZeros[i]...), append(allPoles, gotPoles[i]...)
	}
	assert.Equal(t, sortedRoots([]complex128{-1, 0, 0, 0}), sortedRoots(allZeros))
	assert.Equal(t, sortedRoots(append(poles, 0)), sortedRoots(allPoles))
}

// The first section takes the pole closest to the unit circle, with its conjugate.
func TestPairingStartsNearestTheUnitCircle(t *testing.T) {
	inner := complex(0.3, 0.2)
	outer := cmplx.Rect(0.95, 1)
	poles := []complex128{inner, cmplx.Conj(inner), outer, cmplx.Conj(outer)}
	zeros := []complex128{-1, -1, 1, 1}
	pr, sections, err := NewPairing(zeros, poles)
	require.NoError(t, err)
	_, gotPoles := drain(pr, sections)
	assert.Equal(t, sortedRoots([]complex128{outer, cmplx.Conj(outer)}), sortedRoots(gotPoles[0]))
}

// A complex root without its conjugate can't make a real section.
func TestPairingRejectsUnpairedConjugate(t *testing.T) {
	_, _, err := NewPairing(nil, []complex128{complex(0.5, 0.5)})
	assert.ErrorIs(t, err, ErrUnpairedConjugate)
	_, _, err = NewPairing([]complex128{complex(0, 1)}, []complex128{0.5})
	assert.ErrorIs(t, err, ErrUnpairedConjugate)
}
