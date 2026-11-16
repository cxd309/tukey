package dsp

import (
	"math"
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
		{"unknown padding", b, a, x, []FiltOption{PadType(Padding(7))}, ErrInvalidOption},
		{"negative padding length", b, a, x, []FiltOption{PadLen(-1)}, ErrInvalidParameter},
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

// as for FiltFilt: steady-state starts and odd padding leave a constant untouched
func TestSOSFiltFiltPreservesConstant(t *testing.T) {
	sos := testSOS(t, 3, Lowpass(0.3))
	x := slices.Repeat([]float64{3}, 50)

	y, err := SOSFiltFilt(sos, x)
	require.NoError(t, err)
	reference.AssertClose(t, "y", y, x, reference.Tolerance{Rel: 1e-14, Scaled: 1e-14})
}

// MATLAB's filtfilt counts each polynomial up to its last nonzero coefficient,
// so trailing zeros don't add to the order, and a leading zero, a delay, does.
func TestEffectiveOrder(t *testing.T) {
	cases := []struct {
		name string
		b, a []float64
		want int
	}{
		{"biquad", []float64{1, 2, 1}, []float64{1, 0.5, 0.2}, 2},
		{"trailing zeros in both", []float64{0.3, 0.2, 0}, []float64{1, -0.5, 0}, 1},
		{"b longer than a", []float64{0.25, 0.5, 0.25}, []float64{1, -0.4, 0}, 2},
		{"FIR with a trailing zero", []float64{0.25, 0.5, 0.25, 0}, []float64{1}, 2},
		{"delay", []float64{0, 0, 1}, []float64{1}, 2},
		{"gain", []float64{2}, []float64{1}, 0},
		{"gain section", []float64{2, 0, 0}, []float64{1, 0, 0}, 0},
		{"zero b", []float64{0, 0}, []float64{1}, 0},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			assert.Equal(t, c.want, effectiveOrder(c.b, c.a))
		})
	}
}

// MATLABPadLen is 3 × MATLAB's order, at least 1, so it gives exactly what that
// explicit PadLen does, for (b, a) with trailing zeros and for sections. An FIR
// filter or a gain forgets its start within its order, so its output can't show
// the padding length, but the shortest x it allows, padLen + 1, does.
func TestMATLABPadLenIsThreeTimesTheOrder(t *testing.T) {
	x := testSignal(60)
	filters := []struct {
		name   string
		b, a   []float64
		padLen int
	}{
		{"trailing zeros", []float64{0.3, 0.2, 0}, []float64{1, -0.5, 0}, 3},
		{"FIR with a trailing zero", []float64{0.25, 0.5, 0.25, 0}, []float64{1}, 6},
		{"gain: at least 1", []float64{2}, []float64{1}, 1},
	}
	for _, f := range filters {
		t.Run(f.name, func(t *testing.T) {
			matlab, err := FiltFilt(f.b, f.a, x, MATLABPadLen())
			require.NoError(t, err)
			explicit, err := FiltFilt(f.b, f.a, x, PadLen(f.padLen))
			require.NoError(t, err)
			assert.Equal(t, explicit, matlab)

			_, err = FiltFilt(f.b, f.a, x[:f.padLen], MATLABPadLen())
			assert.ErrorIs(t, err, ErrSignalTooShort, "x as long as the padding")
			_, err = FiltFilt(f.b, f.a, x[:f.padLen+1], MATLABPadLen())
			assert.NoError(t, err, "x one longer than the padding")
		})
	}

	t.Run("odd-order sections", func(t *testing.T) {
		sos := testSOS(t, 3, Lowpass(0.3)) // order 4 to MATLAB: 12, not 9
		matlab, err := SOSFiltFilt(sos, x, BySection(), MATLABPadLen())
		require.NoError(t, err)
		explicit, err := SOSFiltFilt(sos, x, BySection(), PadLen(12))
		require.NoError(t, err)
		assert.Equal(t, explicit, matlab)
	})
}

// BySection is FiltFilt on each section in turn, each with the whole filter's
// padding length, so it's exactly that.
func TestBySectionIsEachSectionInTurn(t *testing.T) {
	sos := testSOS(t, 2, Bandstop(0.2, 0.5))
	x := testSignal(100)
	const padLen = 12 // MATLAB's, for order 4

	want := x
	for _, section := range sos {
		var err error
		want, err = FiltFilt(section[:3], section[3:], want, PadLen(padLen))
		require.NoError(t, err)
	}
	got, err := SOSFiltFilt(sos, x, BySection(), MATLABPadLen())
	require.NoError(t, err)
	assert.Equal(t, want, got)
}

// A (b, a) filter, or one section, is a cascade of one, so filtering it by
// section changes nothing.
func TestBySectionOfOneSectionChangesNothing(t *testing.T) {
	x := testSignal(50)
	sos := testSOS(t, 2, Lowpass(0.3))
	b, a := sos[0][:3], sos[0][3:]

	for _, opts := range [][]FiltOption{nil, {MATLABPadLen()}, {PadType(PaddingEven), PadLen(10)}} {
		plain, err := FiltFilt(b, a, x, opts...)
		require.NoError(t, err)
		bySection, err := FiltFilt(b, a, x, append(slices.Clone(opts), BySection())...)
		require.NoError(t, err)
		assert.Equal(t, plain, bySection)

		plainSOS, err := SOSFiltFilt(sos, x, opts...)
		require.NoError(t, err)
		bySectionSOS, err := SOSFiltFilt(sos, x, append(slices.Clone(opts), BySection())...)
		require.NoError(t, err)
		assert.Equal(t, plainSOS, bySectionSOS)
	}
}

// BySection differs from the cascade only by its edge transients: for a fast
// filter, they agree to rounding away from the ends, and differ at them.
func TestBySectionDiffersFromTheCascadeAtTheEnds(t *testing.T) {
	sos := testSOS(t, 3, Lowpass(0.3))
	x := testSignal(400)

	cascade, err := SOSFiltFilt(sos, x)
	require.NoError(t, err)
	bySection, err := SOSFiltFilt(sos, x, BySection(), MATLABPadLen())
	require.NoError(t, err)

	reference.AssertClose(t, "middle", bySection[100:300], cascade[100:300], reference.Tolerance{Scaled: 1e-14})
	assert.Greater(t, maxAbsDifference(bySection[:3], cascade[:3]), 1e-6, "the ends should differ")
}

// maxAbsDifference is the largest |a[i] − b[i]|
func maxAbsDifference(a, b []float64) (d float64) {
	for i := range a {
		d = max(d, math.Abs(a[i]-b[i]))
	}
	return
}

// Each section starts from its own steady state, so a constant passes a
// unity-DC-gain filter unchanged, ends included.
func TestBySectionPreservesConstant(t *testing.T) {
	sos := testSOS(t, 5, Lowpass(0.2))
	x := slices.Repeat([]float64{3}, 60)
	y, err := SOSFiltFilt(sos, x, BySection(), MATLABPadLen())
	require.NoError(t, err)
	reference.AssertClose(t, "y", y, x, reference.Tolerance{Rel: 1e-14, Scaled: 1e-14})
}

// The options combine in any order.
func TestBySectionOptionOrder(t *testing.T) {
	sos := testSOS(t, 4, Bandpass(0.2, 0.4))
	x := testSignal(80)
	first, err := SOSFiltFilt(sos, x, BySection(), MATLABPadLen(), PadType(PaddingEven))
	require.NoError(t, err)
	last, err := SOSFiltFilt(sos, x, PadType(PaddingEven), MATLABPadLen(), BySection())
	require.NoError(t, err)
	assert.Equal(t, first, last)
}

// The whole filter's padding length decides how short x can be, for every
// section: butter 3 is order 4 to MATLAB, so x needs 13 samples, not 10.
func TestBySectionPaddingLength(t *testing.T) {
	sos := testSOS(t, 3, Lowpass(0.3))
	x := testSignal(20)

	_, err := SOSFiltFilt(sos, x[:12], BySection(), MATLABPadLen())
	assert.ErrorIs(t, err, ErrSignalTooShort)
	_, err = SOSFiltFilt(sos, x[:13], BySection(), MATLABPadLen())
	assert.NoError(t, err)

	// as without BySection, PaddingNone ignores the padding length
	_, err = SOSFiltFilt(sos, x[:3], BySection(), PadType(PaddingNone), PadLen(100))
	assert.NoError(t, err)
}

func TestBySectionRejectsInvalidInput(t *testing.T) {
	sos := testSOS(t, 2, Lowpass(0.3))
	x := testSignal(40)
	cases := []struct {
		name    string
		sos     SOS
		opts    []FiltOption
		wantErr error
	}{
		{"unknown padding", sos, []FiltOption{PadType(Padding(7))}, ErrInvalidOption},
		{"negative padding length", sos, []FiltOption{PadLen(-1)}, ErrInvalidParameter},
		{"pole at DC", SOS{{1, 0, 0, 1, -1, 0}}, nil, ErrNoSteadyState},
		{"pole at DC in a later section", append(slices.Clone(sos), [6]float64{1, 0, 0, 1, -1, 0}), nil, ErrNoSteadyState},
		{"invalid sections", SOS{{1, 0, 0, 0, 0, 0}}, nil, ErrInvalidSOS},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			_, err := SOSFiltFilt(c.sos, x, append(slices.Clone(c.opts), BySection())...)
			assert.ErrorIs(t, err, c.wantErr)
		})
	}
}

func TestBySectionDoesNotModifyInputs(t *testing.T) {
	sos := testSOS(t, 3, Lowpass(0.3))
	x := testSignal(40)
	sosBefore, xBefore := slices.Clone(sos), slices.Clone(x)
	_, err := SOSFiltFilt(sos, x, BySection(), MATLABPadLen())
	require.NoError(t, err)
	assert.Equal(t, sosBefore, sos)
	assert.Equal(t, xBefore, x)
}

func TestPaddingExtend(t *testing.T) {
	x := []float64{1, 2, 4, 7}
	cases := []struct {
		name    string
		padding Padding
		n       int
		want    []float64
	}{
		// left: 2*1-4, 2*1-2   right: 2*7-4, 2*7-2
		{"odd", PaddingOdd, 2, []float64{-2, 0, 1, 2, 4, 7, 10, 12}},
		{"even", PaddingEven, 2, []float64{4, 2, 1, 2, 4, 7, 4, 2}},
		{"constant", PaddingConstant, 2, []float64{1, 1, 1, 2, 4, 7, 7, 7}},
		{"longest odd extension, n = len(x)-1", PaddingOdd, 3, []float64{-5, -2, 0, 1, 2, 4, 7, 10, 12, 13}},
		{"none", PaddingNone, 0, []float64{1, 2, 4, 7}},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			xBefore := slices.Clone(x)
			assert.Equal(t, c.want, c.padding.extend(x, c.n))
			assert.Equal(t, xBefore, x, "extend must not modify x")
		})
	}
}

func TestPaddingExtendReturnsNewSlice(t *testing.T) {
	x := []float64{1, 2, 3}
	extended := PaddingNone.extend(x, 0)
	extended[0] = 99
	assert.Equal(t, 1.0, x[0], "extend must not alias x")
}

func TestPaddingExtendPanicsOnBrokenPreconditions(t *testing.T) {
	cases := []struct {
		name    string
		padding Padding
		x       []float64
		n       int
	}{
		{"n = len(x)", PaddingOdd, []float64{1, 2, 3}, 3},
		{"negative n", PaddingOdd, []float64{1, 2, 3}, -1},
		{"empty x", PaddingOdd, nil, 0},
		{"PaddingNone with n > 0", PaddingNone, []float64{1, 2, 3}, 1},
		{"unknown padding", Padding(7), []float64{1, 2, 3}, 1},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			assert.Panics(t, func() { c.padding.extend(c.x, c.n) })
		})
	}
}
