package window

import (
	"math"
	"testing"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

// GeneralGaussian at p = 1 is Gaussian, to within rounding: one squares, the other uses Pow
func TestGeneralGaussianAtP1IsGaussian(t *testing.T) {
	for n := 0; n <= 20; n++ {
		gaussian, err := Gaussian(n, 2.5)
		require.NoError(t, err)
		general, err := GeneralGaussian(n, 1, 2.5)
		require.NoError(t, err)
		assert.InDeltaSlice(t, gaussian, general, 1e-15, "n = %d", n)
	}
}

// Center at the periodic window's own middle (n/2) is the default, exactly:
// |i - n/2| is the same whether computed directly or mirrored
func TestExponentialCenterAtMiddleIsDefault(t *testing.T) {
	for n := 0; n <= 20; n++ {
		centred, err := Exponential(n, 3, Periodic())
		require.NoError(t, err)
		w, err := Exponential(n, 3, Periodic(), Center(float64(n)/2))
		require.NoError(t, err)
		assert.Equal(t, centred, w, "n = %d", n)
	}
}

func TestExponentialFamilyRejectsInvalidParameters(t *testing.T) {
	assertInvalidParameter(t, map[string]func() ([]float64, error){
		"Gaussian std 0":        func() ([]float64, error) { return Gaussian(8, 0) },
		"Gaussian std NaN":      func() ([]float64, error) { return Gaussian(8, math.NaN()) },
		"GeneralGaussian p -1":  func() ([]float64, error) { return GeneralGaussian(8, -1, 2) },
		"GeneralGaussian sig 0": func() ([]float64, error) { return GeneralGaussian(8, 1, 0) },
		"Exponential tau 0":     func() ([]float64, error) { return Exponential(8, 0) },
		"Exponential Center without Periodic": func() ([]float64, error) {
			return Exponential(8, 3, Center(0))
		},
		"Exponential Center NaN": func() ([]float64, error) {
			return Exponential(8, 3, Periodic(), Center(math.NaN()))
		},
	})
}
