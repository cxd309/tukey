package window

import (
	"math"
	"testing"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

// Cosine is half a period of a sine, sampled at the middle of each of n equal steps,
// so it never reaches 0 at its ends, which are sin(π/(2n)).
func TestCosineIsHalfASine(t *testing.T) {
	for _, n := range []int{2, 7, 8, 64} {
		w, err := Cosine(n)
		require.NoError(t, err)
		for k, v := range w {
			assert.InDelta(t, math.Sin(math.Pi*(float64(k)+0.5)/float64(n)), v, 1e-15, "n %d, k %d", n, k)
		}
		assert.Greater(t, w[0], 0.0, "n %d", n)
	}
}

// Lanczos is the central lobe of a sinc, sin(π·x)/(π·x) for x from -1 to 1,
// so it's 1 at the centre of an odd n and 0 at both ends, rounding apart.
func TestLanczosIsTheSincsCentralLobe(t *testing.T) {
	for _, n := range []int{3, 8, 9, 65} {
		w, err := Lanczos(n)
		require.NoError(t, err)
		for k, v := range w {
			x := 2*float64(k)/float64(n-1) - 1
			want := 1.0
			if x != 0 {
				want = math.Sin(math.Pi*x) / (math.Pi * x)
			}
			assert.InDelta(t, want, v, 1e-15, "n %d, k %d", n, k)
		}
		assert.InDelta(t, 0, w[0], 1e-15, "n %d", n)
		if n%2 == 1 {
			assert.Equal(t, 1.0, w[n/2], "n %d", n)
		}
	}
}
