package fft

import (
	"strconv"
	"testing"

	"github.com/cxd309/tukey/internal/reference"
	"github.com/stretchr/testify/assert"
)

// Each butterfly is the p-point DFT:
// the dedicated ones for 2, 3, 4 and 5, the direct sum for the primes up to maxRadix, and Bluestein above it.
// Each writes its outputs step apart, so they're read back from every other slot.
func TestButterflyIsDFT(t *testing.T) {
	const step = 2
	for _, p := range []int{2, 3, 4, 5, 7, 11, 13, 17} {
		t.Run(strconv.Itoa(p), func(t *testing.T) {
			terms := randomSignal(p)
			dst := make([]complex128, p*step)
			butterfly(terms, twiddleTable(p), 1, dst, step)
			got := make([]complex128, p)
			for s := range got {
				got[s] = dst[s*step]
			}
			reference.AssertClose(t, "X", got, dft(terms), reference.Tolerance{Rel: 1e-14, Scaled: 1e-14})
		})
	}
}

// minusI rotates by a quarter turn exactly, and scaledBy multiplies each part once.
func TestRotationsAreExact(t *testing.T) {
	re, im, k := 0.1, -3.7, 0.3
	z := complex(re, im)
	assert.Equal(t, complex(im, -re), minusI(z))
	assert.Equal(t, complex(k*re, k*im), scaledBy(k, z))
}
