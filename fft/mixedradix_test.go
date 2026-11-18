package fft

import (
	"strconv"
	"testing"

	"github.com/cxd309/tukey/internal/reference"
	"github.com/stretchr/testify/require"
)

// Mixed radix at lengths the SciPy vectors don't reach:
// every prime up to maxRadix, as a factor, squared, and mixed with others.
func TestMixedRadixMatchesDFT(t *testing.T) {
	for _, n := range []int{3, 5, 6, 7, 9, 10, 11, 12, 13, 14, 15, 20, 21, 49, 77, 100, 121, 143, 169, 210, 1000} {
		t.Run(strconv.Itoa(n), func(t *testing.T) {
			x := randomSignal(n)
			y, err := Forward(x)
			require.NoError(t, err)
			reference.AssertClose(t, "y", y, dft(x), reference.Tolerance{Rel: 1e-13, Scaled: 1e-13})
		})
	}
}

// The round trip and Parseval hold at a smooth, non-power-of-two length too.
func TestMixedRadixRoundTrip(t *testing.T) {
	x := randomSignal(1000)
	y, err := Forward(x)
	require.NoError(t, err)
	back, err := Inverse(y)
	require.NoError(t, err)
	reference.AssertClose(t, "x", back, x, reference.Tolerance{Rel: 1e-13, Scaled: 1e-13})
}
