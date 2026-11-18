package fft

import (
	"strconv"
	"testing"

	"github.com/cxd309/tukey/internal/reference"
	"github.com/stretchr/testify/require"
)

// The fast transform must agree with the definition, at every power of two.
func TestRadix2MatchesDFT(t *testing.T) {
	for n := 1; n <= 256; n *= 2 {
		t.Run(strconv.Itoa(n), func(t *testing.T) {
			x := randomSignal(n)
			y, err := Forward(x)
			require.NoError(t, err)
			reference.AssertClose(t, "y", y, dft(x), reference.Tolerance{Rel: 1e-13, Scaled: 1e-13})
		})
	}
}

// The inverse undoes the forward transform.
func TestRadix2RoundTrip(t *testing.T) {
	x := randomSignal(1024)
	y, err := Forward(x)
	require.NoError(t, err)
	back, err := Inverse(y)
	require.NoError(t, err)
	reference.AssertClose(t, "x", back, x, reference.Tolerance{Rel: 1e-13, Scaled: 1e-13})
}
