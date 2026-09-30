package dsp

import (
	"slices"
	"testing"

	"github.com/cxd309/tukey/internal/reference"
	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

func TestFilterRejectsInvalidInput(t *testing.T) {
	cases := []struct {
		name    string
		b, a    []float64
		wantErr error
	}{
		{"empty b", nil, []float64{1}, ErrInvalidCoefficients},
		{"empty a", []float64{1}, nil, ErrInvalidCoefficients},
		{"a[0] zero", []float64{1}, []float64{0, 1}, ErrInvalidCoefficients},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			_, err := Filter(c.b, c.a, []float64{1, 2, 3})
			assert.ErrorIs(t, err, c.wantErr)
		})
	}
}

// test never mutate inputs
// unnormalised a[0] != 1 case is the one where normalising in place would be tempting
func TestFilterDoesNotModifyInputs(t *testing.T) {
	b, a, x := []float64{2, 1}, []float64{2, -0.5}, []float64{1, 2, 3}
	bBefore, aBefore, xBefore := slices.Clone(b), slices.Clone(a), slices.Clone(x)

	_, err := Filter(b, a, x)
	require.NoError(t, err)

	assert.Equal(t, bBefore, b)
	assert.Equal(t, aBefore, a)
	assert.Equal(t, xBefore, x)
}

// the point of steadyState: a filter started there treats a constant input as
// already settled, so the output is the input times the DC gain from the very first
// sample, with none of the ramp and overshoot seen when starting from rest
func TestSteadyStateRemovesStartupTransient(t *testing.T) {
	fButter, err := Butter(2, Lowpass(0.3)) // unity DC gain
	require.NoError(t, err)
	b, a, err := fButter.BA()
	require.NoError(t, err)
	f, err := newDigitalFilter(b, a)
	require.NoError(t, err)
	zi, err := f.steadyState()
	require.NoError(t, err)

	const level = 3.0
	for i := range zi {
		zi[i] *= level
	}
	x := slices.Repeat([]float64{level}, 50)

	y, _ := f.apply(x, zi)
	reference.AssertClose(t, "y", y, x, reference.Tolerance{Rel: 1e-14, Scaled: 1e-14})
}

func TestSteadyStateRejectsPoleAtDC(t *testing.T) {
	f, err := newDigitalFilter([]float64{1}, []float64{1, -1}) // integrator
	require.NoError(t, err)
	_, err = f.steadyState()
	assert.ErrorIs(t, err, ErrNoSteadyState)
}

// as for SOSFilterState: two blocks with carried state must match one call exactly
func TestFilterStateBlocksMatchWhole(t *testing.T) {
	f, err := Butter(4, Bandpass(0.2, 0.5))
	require.NoError(t, err)
	b, a, err := f.BA()
	require.NoError(t, err)
	x := testSignal(200)

	whole, err := Filter(b, a, x)
	require.NoError(t, err)
	first, zf, err := FilterState(b, a, x[:73], make([]float64, len(a)-1))
	require.NoError(t, err)
	second, _, err := FilterState(b, a, x[73:], zf)
	require.NoError(t, err)

	assert.Equal(t, whole, slices.Concat(first, second))
}

func TestFilterStateRejectsWrongStateLength(t *testing.T) {
	_, _, err := FilterState([]float64{1, 1}, []float64{1, -0.5}, []float64{1, 2, 3}, []float64{0, 0})
	assert.ErrorIs(t, err, ErrInvalidState, "two taps need one state value")
}
