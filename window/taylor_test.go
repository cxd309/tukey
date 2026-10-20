package window

import (
	"testing"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

// with nbar = 1 there are no cosine terms: Taylor is Boxcar, exactly
func TestTaylorAtNbar1IsBoxcar(t *testing.T) {
	for n := 0; n <= 20; n++ {
		boxcar, err := Boxcar(n)
		require.NoError(t, err)
		w, err := Taylor(n, 1, 30)
		require.NoError(t, err)
		assert.Equal(t, boxcar, w, "n = %d", n)
	}
}

// normalising divides every sample by the centre's value, so for odd n the
// normalised window is the unnormalised one over its own centre, exactly
func TestTaylorNormalisesToCentre(t *testing.T) {
	for _, n := range []int{3, 9, 65} {
		normalised, err := Taylor(n, 4, 30)
		require.NoError(t, err)
		raw, err := Taylor(n, 4, 30, Unnormalised())
		require.NoError(t, err)
		centre := raw[n/2]
		for k := range raw {
			assert.Equal(t, raw[k]/centre, normalised[k], "n = %d, k = %d", n, k)
		}
		assert.Equal(t, 1.0, normalised[n/2])
	}
}

func TestTaylorRejectsInvalidParameters(t *testing.T) {
	assertInvalidParameter(t, map[string]func() ([]float64, error){
		"nbar 0": func() ([]float64, error) { return Taylor(8, 0, 30) },
		"sll 0":  func() ([]float64, error) { return Taylor(8, 4, 0) },
	})
}
