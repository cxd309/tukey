package fft

import (
	"strconv"
	"testing"

	"github.com/cxd309/tukey/internal/reference"
	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

// Bluestein, on its own for large primes, and as the butterfly inside mixed radix.
func TestBluesteinMatchesDFT(t *testing.T) {
	for _, n := range []int{17, 19, 23, 97, 101, 34, 291, 1009, 2018} {
		t.Run(strconv.Itoa(n), func(t *testing.T) {
			x := randomSignal(n)
			y, err := Forward(x)
			require.NoError(t, err)
			reference.AssertClose(t, "y", y, dft(x), reference.Tolerance{Rel: 1e-13, Scaled: 1e-13})
		})
	}
}

// The round trip holds through Bluestein, at a prime length.
func TestBluesteinRoundTrip(t *testing.T) {
	x := randomSignal(1009)
	y, err := Forward(x)
	require.NoError(t, err)
	back, err := Inverse(y)
	require.NoError(t, err)
	reference.AssertClose(t, "x", back, x, reference.Tolerance{Rel: 1e-13, Scaled: 1e-13})
}

func TestSmoothSize(t *testing.T) {
	cases := map[int]int{
		1: 1, 7: 8, 11: 12, 13: 15, 17: 18,
		2017:   2025,   // 3^4 * 5^2, not 2048
		131073: 131220, // 2^2 * 3^8 * 5, half of 262144
	}
	for target, want := range cases {
		assert.Equal(t, want, smoothSize(target), "smoothSize(%d)", target)
	}
}

// smoothSize must give the smallest smooth length, so check it against a brute-force search.
func TestSmoothSizeIsSmallest(t *testing.T) {
	isSmooth := func(n int) bool {
		for _, f := range []int{2, 3, 5} {
			for n%f == 0 {
				n /= f
			}
		}
		return n == 1
	}
	for target := 1; target <= 2000; target++ {
		want := target
		for !isSmooth(want) {
			want++
		}
		assert.Equal(t, want, smoothSize(target), "smoothSize(%d)", target)
	}
}

// The padded length is the power of two unless the smooth length is enough shorter.
func TestPaddedSize(t *testing.T) {
	assert.Equal(t, 2048, paddedSize(1009))    // 2025 is barely shorter
	assert.Equal(t, 131220, paddedSize(65537)) // 131220 is half of 262144
}
