package dsp

import (
	"slices"
	"testing"

	"github.com/cxd309/tukey/internal/reference"
	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

type lfilterZiVector struct {
	reference.Meta
	Params struct {
		B []float64 `json:"b"`
		A []float64 `json:"a"`
	} `json:"params"`
	Output struct {
		Zi []float64 `json:"zi"`
	} `json:"output"`
}

// steadyStateTolerance is a placeholder: measure it, then replace this comment
var steadyStateTolerance = reference.Tolerance{Rel: 1e-14, Scaled: 1e-14}

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

func TestSteadyStateMatchesLfilterZi(t *testing.T) {
	reference.Run(t, "lfilter_zi", func(t *testing.T, v lfilterZiVector) {
		f, err := newDigitalFilter(v.Params.B, v.Params.A)
		require.NoError(t, err)
		zi, err := f.steadyState()
		require.NoError(t, err)
		reference.AssertClose(t, "zi", zi, v.Output.Zi, steadyStateTolerance)
	})
}

// the point of steadyState: a filter started there treats a constant input as
// already settled, so the output is the input times the DC gain from the very first
// sample, with none of the ramp and overshoot seen when starting from rest
func TestSteadyStateRemovesStartupTransient(t *testing.T) {
	b, a, err := Butter(2, Lowpass(0.3)) // unity DC gain
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

	y := f.apply(x, zi)
	reference.AssertClose(t, "y", y, x, reference.Tolerance{Rel: 1e-14, Scaled: 1e-14})
}

func TestSteadyStateRejectsPoleAtDC(t *testing.T) {
	f, err := newDigitalFilter([]float64{1}, []float64{1, -1}) // integrator
	require.NoError(t, err)
	_, err = f.steadyState()
	assert.ErrorIs(t, err, ErrNoSteadyState)
}
