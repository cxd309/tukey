package dsp

import (
	"slices"
	"testing"

	"github.com/cxd309/tukey/internal/reference"
	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

func TestFiltFiltRejectsInvalidInput(t *testing.T) {
	b, a := []float64{0.5, 0.5}, []float64{1, -0.5} // 2 taps: default padding length 6
	x := slices.Repeat([]float64{1}, 20)
	cases := []struct {
		name    string
		b, a, x []float64
		opts    []FiltOption
		wantErr error
	}{
		{"empty b", nil, a, x, nil, ErrInvalidCoefficients},
		{"unknown padding", b, a, x, []FiltOption{PadType(Padding(7))}, ErrInvalidPadding},
		{"negative padding length", b, a, x, []FiltOption{PadLen(-1)}, ErrInvalidPadding},
		{"x as long as the padding", b, a, x[:6], nil, ErrSignalTooShort},
		{"empty x", b, a, nil, []FiltOption{PadType(PaddingNone)}, ErrSignalTooShort},
		{"pole at DC", b, []float64{1, -1}, x, nil, ErrNoSteadyState},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			_, err := FiltFilt(c.b, c.a, c.x, c.opts...)
			assert.ErrorIs(t, err, c.wantErr)
		})
	}
}

func TestFiltFiltDoesNotModifyInputs(t *testing.T) {
	b, a := []float64{2, 1}, []float64{2, -0.5}
	x := []float64{1, 2, 4, 7, 11, 16, 22, 29}
	bBefore, aBefore, xBefore := slices.Clone(b), slices.Clone(a), slices.Clone(x)

	_, err := FiltFilt(b, a, x)
	require.NoError(t, err)

	assert.Equal(t, bBefore, b)
	assert.Equal(t, aBefore, a)
	assert.Equal(t, xBefore, x)
}

// with steady-state initial conditions and odd padding, a constant input has no
// transient anywhere, so a unity-DC-gain filter returns it unchanged at every
// sample, ends included, to within rounding (SciPy's own output is 1 ULP off)
func TestFiltFiltPreservesConstant(t *testing.T) {
	f, err := Butter(2, Lowpass(0.3))
	require.NoError(t, err)
	b, a, err := f.BA()
	require.NoError(t, err)
	x := slices.Repeat([]float64{3}, 50)

	y, err := FiltFilt(b, a, x)
	require.NoError(t, err)
	reference.AssertClose(t, "y", y, x, reference.Tolerance{Rel: 1e-14, Scaled: 1e-14})
}

// MATLABPadLen is only a padding length, so it must give exactly
// what the equivalent explicit PadLen does: 3 * (3 taps - 1) = 6
func TestMATLABPadLenMatchesExplicitPadLen(t *testing.T) {
	f, err := Butter(2, Lowpass(0.3))
	require.NoError(t, err)
	b, a, err := f.BA()
	require.NoError(t, err)
	x := []float64{0, 1, 3, 2, 5, 4, 6, 8, 7, 9, 10, 12}

	matlab, err := FiltFilt(b, a, x, MATLABPadLen())
	require.NoError(t, err)
	explicit, err := FiltFilt(b, a, x, PadLen(6))
	require.NoError(t, err)
	assert.Equal(t, explicit, matlab)
}

// as in SciPy, PaddingNone means nothing is added, so PadLen is ignored
// and x only has to be non-empty
func TestFiltFiltPaddingNoneIgnoresPadLen(t *testing.T) {
	b, a := []float64{0.5, 0.5}, []float64{1, -0.5}
	_, err := FiltFilt(b, a, []float64{1, 2, 3}, PadType(PaddingNone), PadLen(100))
	assert.NoError(t, err)
}
