package dsp

import (
	"math"
	"math/cmplx"
	"testing"

	"github.com/cxd309/tukey/internal/reference"
	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

// the defining property of a Butterworth lowpass: gain 1 at DC, and exactly
// 1/sqrt(2) (-3 dB) at the cutoff, which prewarping puts precisely where asked;
// with n = 10, w[3] is 0.3*pi, the cutoff
func TestFreqzButterworthIsMinus3dBAtCutoff(t *testing.T) {
	f, err := Butter(4, Lowpass(0.3))
	require.NoError(t, err)
	b, a, err := f.BA()
	require.NoError(t, err)

	w, h, err := Freqz(b, a, 10)
	require.NoError(t, err)
	assert.InDelta(t, 0.3*math.Pi, w[3], 1e-15)
	assert.InDelta(t, 1, cmplx.Abs(h[0]), 1e-12, "DC gain")
	assert.InDelta(t, 1/math.Sqrt2, cmplx.Abs(h[3]), 1e-12, "gain at the cutoff")
}

// sections and (b, a) are the same filter, so at a well-conditioned order
// they must have the same response
func TestSOSFreqzMatchesFreqz(t *testing.T) {
	f, err := Butter(4, Bandpass(0.2, 0.5))
	require.NoError(t, err)
	b, a, err := f.BA()
	require.NoError(t, err)
	sos, err := f.SOS()
	require.NoError(t, err)

	_, viaBA, err := Freqz(b, a, 64)
	require.NoError(t, err)
	_, viaSOS, err := SOSFreqz(sos, 64)
	require.NoError(t, err)
	reference.AssertClose(t, "h", viaSOS, viaBA, reference.Tolerance{Rel: 1e-12, Scaled: 1e-12})
}

func TestFreqzRejectsInvalidInput(t *testing.T) {
	_, _, err := Freqz(nil, []float64{1}, 8)
	assert.ErrorIs(t, err, ErrInvalidCoefficients)
	_, _, err = Freqz([]float64{1}, []float64{1}, -1)
	assert.ErrorIs(t, err, ErrInvalidLength)
	_, _, err = SOSFreqz(SOS{{1, 0, 0, 2, 0, 0}}, 8)
	assert.ErrorIs(t, err, ErrInvalidSOS)
}
