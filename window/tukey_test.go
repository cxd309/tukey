package window

import (
	"math"
	"testing"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

// Tukey's alpha has no domain: at or below 0 it's Boxcar, at or above 1 it's Hann
func TestTukeyLimitsAreExact(t *testing.T) {
	for n := 0; n <= 20; n++ {
		boxcar, err := Boxcar(n)
		require.NoError(t, err)
		hann, err := Hann(n)
		require.NoError(t, err)

		for _, alpha := range []float64{-1, 0} {
			w, err := Tukey(n, alpha)
			require.NoError(t, err)
			assert.Equal(t, boxcar, w, "Tukey(%d, %g)", n, alpha)
		}
		for _, alpha := range []float64{1, 2} {
			w, err := Tukey(n, alpha)
			require.NoError(t, err)
			assert.Equal(t, hann, w, "Tukey(%d, %g)", n, alpha)
		}
	}
}

func TestTukeyRejectsNaN(t *testing.T) {
	_, err := Tukey(8, math.NaN())
	assert.ErrorIs(t, err, ErrInvalidParameter)
}
