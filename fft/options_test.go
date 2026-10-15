package fft

import (
	"slices"
	"testing"

	"github.com/cxd309/tukey/internal/reference"
	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

// whichever normalisation, a forward and inverse pair recovers the signal
func TestEveryNormalisationRoundTrips(t *testing.T) {
	x := randomSignal(100)
	for _, norm := range []Normalisation{NormBackward, NormOrtho, NormForward} {
		y, err := Forward(x, Norm(norm))
		require.NoError(t, err)
		back, err := Inverse(y, Norm(norm))
		require.NoError(t, err)
		reference.AssertClose(t, "x", back, x, reference.Tolerance{Rel: 1e-13, Scaled: 1e-13})
	}
}

// NormOrtho makes the transform unitary: energy is identical in both domains,
// with no 1/n correction (compare TestForwardPreservesEnergy)
func TestOrthoPreservesEnergyDirectly(t *testing.T) {
	x := randomSignal(256)
	y, err := Forward(x, Norm(NormOrtho))
	require.NoError(t, err)
	timeEnergy, freqEnergy := 0.0, 0.0
	for i := range x {
		timeEnergy += real(x[i])*real(x[i]) + imag(x[i])*imag(x[i])
		freqEnergy += real(y[i])*real(y[i]) + imag(y[i])*imag(y[i])
	}
	assert.InEpsilon(t, timeEnergy, freqEnergy, 1e-13)
}

// the default is NormBackward, bit for bit
func TestDefaultNormalisationIsBackward(t *testing.T) {
	x := randomSignal(60)
	plain, err := Forward(x)
	require.NoError(t, err)
	explicit, err := Forward(x, Norm(NormBackward))
	require.NoError(t, err)
	assert.Equal(t, plain, explicit)
}

// N zero-pads exactly as padding by hand would
func TestNPadsWithZeros(t *testing.T) {
	x := randomSignal(12)
	viaN, err := Forward(x, N(16))
	require.NoError(t, err)
	byHand, err := Forward(append(slices.Clone(x), make([]complex128, 4)...))
	require.NoError(t, err)
	assert.Equal(t, byHand, viaN)
}

// shifting the frequencies puts them in ascending order, and InverseShift undoes Shift
func TestShiftOrdersFrequencies(t *testing.T) {
	for _, n := range []int{8, 9} {
		f := Freq(n, 1)
		assert.True(t, slices.IsSorted(Shift(f)), "n=%d", n)
		assert.Equal(t, f, InverseShift(Shift(f)), "n=%d", n)
	}
}

func TestOptionsRejectInvalidInput(t *testing.T) {
	x := randomSignal(8)
	_, err := Forward(x, Norm(Normalisation(7)))
	assert.ErrorIs(t, err, ErrInvalidOption, "unknown normalisation")
	_, err = Forward(x, N(0))
	assert.ErrorIs(t, err, ErrInvalidLength, "N(0)")
	_, err = InverseReal([]complex128{1, 2}, 2, N(4))
	assert.ErrorIs(t, err, ErrInvalidOption, "N on InverseReal")
	_, err = Forward(nil)
	assert.ErrorIs(t, err, ErrInvalidLength, "empty input")
}
