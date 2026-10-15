package fft

import (
	"testing"

	"github.com/stretchr/testify/assert"
)

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

// smoothSize must give the smallest smooth length, so check it against a brute-force search
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

// the padded length is the power of two unless the smooth length is enough shorter
func TestPaddedSize(t *testing.T) {
	assert.Equal(t, 2048, paddedSize(1009))    // 2025 is barely shorter
	assert.Equal(t, 131220, paddedSize(65537)) // 131220 is half of 262144
}
