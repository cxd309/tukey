package window

import (
	"math"
	"math/cmplx"
	"testing"

	"github.com/cxd309/tukey/fft"
	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

// the property Chebwin is designed for: every sidelobe of its spectrum is at
// -attenuation dB relative to the main lobe; checked on a finely sampled spectrum
func TestChebwinSidelobesAreAtAttenuation(t *testing.T) {
	for _, attenuation := range []float64{50, 80} {
		w, err := Chebwin(31, attenuation)
		require.NoError(t, err)
		spectrum, err := fft.ForwardReal(w, fft.N(8192))
		require.NoError(t, err)

		// walk down the main lobe to its first null; everything after is sidelobes
		k := 1
		for k < len(spectrum)-1 && cmplx.Abs(spectrum[k+1]) < cmplx.Abs(spectrum[k]) {
			k++
		}
		peak := 0.0
		for _, s := range spectrum[k:] {
			peak = math.Max(peak, cmplx.Abs(s))
		}
		level := peak / cmplx.Abs(spectrum[0])
		assert.InEpsilon(t, math.Pow(10, -attenuation/20), level, 1e-2, "%g dB", attenuation)
	}
}

func TestChebwinRejectsInvalidParameters(t *testing.T) {
	assertInvalidParameter(t, map[string]func() ([]float64, error){
		"attenuation 0":   func() ([]float64, error) { return Chebwin(8, 0) },
		"attenuation NaN": func() ([]float64, error) { return Chebwin(8, math.NaN()) },
	})
}
