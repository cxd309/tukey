package fft

import (
	"math/cmplx"
	"testing"

	"github.com/stretchr/testify/assert"
)

// quarter and half turns are exact, not merely close: exp(-2πi·k/n) computed
// directly gives e.g. 6.1e-17 - 1i for a quarter turn, where this gives exactly -i
func TestTwiddleExactAtQuarterTurns(t *testing.T) {
	const n = 12
	assert.Equal(t, complex(1, 0), twiddle(0, n))
	assert.Equal(t, complex(0, -1), twiddle(n/4, n))
	assert.Equal(t, complex(-1, 0), twiddle(n/2, n))
	assert.Equal(t, complex(0, 1), twiddle(3*n/4, n))
}

// going round the other way is the conjugate, exactly, because both reduce to
// the same small angle
func TestTwiddleConjugateSymmetry(t *testing.T) {
	const n = 97
	for k := range n {
		assert.Equal(t, cmplx.Conj(twiddle(k, n)), twiddle(n-k, n), "k=%d", k)
	}
}

// k can be negative or beyond n: only k mod n matters
func TestTwiddleReducesK(t *testing.T) {
	assert.Equal(t, twiddle(3, 10), twiddle(13, 10))
	assert.Equal(t, twiddle(7, 10), twiddle(-3, 10))
}
