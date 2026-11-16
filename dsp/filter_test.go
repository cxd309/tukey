package dsp

import (
	"math"
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
		{"NaN in b", []float64{1, math.NaN()}, []float64{1}, ErrInvalidCoefficients},
		{"infinite a[0]", []float64{1}, []float64{math.Inf(1), 1}, ErrInvalidCoefficients},
		{"infinite a[1]", []float64{1}, []float64{1, math.Inf(-1)}, ErrInvalidCoefficients},
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

// testSOS designs a filter as sections for these tests
func testSOS(t *testing.T, order int, band Band) (sos SOS) {
	t.Helper()
	f, err := Butter(order, band)
	require.NoError(t, err)
	sos, err = f.SOS()
	require.NoError(t, err)
	return
}

// testSignal is a deterministic signal with some structure: two sines and a step
func testSignal(n int) (x []float64) {
	x = make([]float64, n)
	for i := range x {
		t := float64(i) / float64(n)
		x[i] = math.Sin(2*math.Pi*5*t) + 0.5*math.Sin(2*math.Pi*40*t)
		if i >= n/2 {
			x[i] += 1
		}
	}
	return
}

// filtering in blocks, carrying each block's final state into the next, must give
// exactly the same output as filtering in one go: it's the same arithmetic in the
// same order, so this is bit-for-bit, not approximate
func TestSOSFilterStateBlocksMatchWhole(t *testing.T) {
	sos := testSOS(t, 5, Lowpass(0.2))
	x := testSignal(200)

	whole, err := SOSFilter(sos, x)
	require.NoError(t, err)

	first, zf, err := SOSFilterState(sos, x[:73], make([][2]float64, len(sos)))
	require.NoError(t, err)
	second, _, err := SOSFilterState(sos, x[73:], zf)
	require.NoError(t, err)

	assert.Equal(t, whole, slices.Concat(first, second))
}

// sections and (b, a) describe the same filter, so at an order where (b, a) is
// still well conditioned they must agree closely
func TestSOSFilterMatchesFilter(t *testing.T) {
	f, err := Butter(4, Bandpass(0.2, 0.5))
	require.NoError(t, err)
	b, a, err := f.BA()
	require.NoError(t, err)
	sos, err := f.SOS()
	require.NoError(t, err)
	x := testSignal(200)

	viaBA, err := Filter(b, a, x)
	require.NoError(t, err)
	viaSOS, err := SOSFilter(sos, x)
	require.NoError(t, err)
	reference.AssertClose(t, "y", viaSOS, viaBA, reference.Tolerance{Rel: 1e-12, Scaled: 1e-12})
}

func TestSOSFilteringRejectsInvalidInput(t *testing.T) {
	valid := SOS{{1, 0, 0, 1, -0.5, 0}}
	x := testSignal(50)

	_, err := SOSFilter(SOS{}, x)
	assert.ErrorIs(t, err, ErrInvalidSOS, "no sections")

	_, err = SOSFilter(SOS{{1, 0, 0, 2, -0.5, 0}}, x)
	assert.ErrorIs(t, err, ErrInvalidSOS, "a0 not 1")

	_, err = SOSFilter(SOS{{1, math.NaN(), 0, 1, -0.5, 0}}, x)
	assert.ErrorIs(t, err, ErrInvalidSOS, "NaN coefficient")

	_, _, err = SOSFilterState(valid, x, [][2]float64{{0, math.Inf(1)}})
	assert.ErrorIs(t, err, ErrInvalidState, "infinite state")

	_, _, err = SOSFilterState(valid, x, make([][2]float64, 2))
	assert.ErrorIs(t, err, ErrInvalidState, "one section, two state pairs")

	_, err = SOSFilterZi(SOS{{1, 0, 0, 1, -1, 0}})
	assert.ErrorIs(t, err, ErrNoSteadyState, "pole at DC")

	_, err = SOSFiltFilt(valid, x[:3])
	assert.ErrorIs(t, err, ErrSignalTooShort, "3 samples, padding length 3*2")
}

func TestSOSFilteringDoesNotModifyInputs(t *testing.T) {
	sos := testSOS(t, 3, Lowpass(0.3))
	x := testSignal(50)
	zi := [][2]float64{{0.1, 0.2}, {0.3, 0.4}}
	sosBefore, xBefore, ziBefore := slices.Clone(sos), slices.Clone(x), slices.Clone(zi)

	_, _, err := SOSFilterState(sos, x, zi)
	require.NoError(t, err)
	_, err = SOSFiltFilt(sos, x)
	require.NoError(t, err)

	assert.Equal(t, sosBefore, sos)
	assert.Equal(t, xBefore, x)
	assert.Equal(t, ziBefore, zi)
}

func TestFilterStateRejectsNonFiniteState(t *testing.T) {
	for _, v := range []float64{math.NaN(), math.Inf(1)} {
		_, _, err := FilterState([]float64{1, 1}, []float64{1, -0.5}, testSignal(10), []float64{v})
		assert.ErrorIs(t, err, ErrInvalidState, "zi %g", v)
	}
}

// A NaN or infinite sample reaches only the outputs whose taps are non-zero,
// where 0·NaN in the state would carry NaN to every later output.
func TestFilterZeroCoefficientsDropNonFinite(t *testing.T) {
	const bad = 5
	cases := []struct {
		name     string
		b        []float64
		affected []int
	}{
		{"FIR", []float64{0.2, 0.4, 0.2}, []int{5, 6, 7}},
		{"FIR with a zero tap", []float64{1, 0, 1}, []int{5, 7}},
	}
	for _, c := range cases {
		for _, v := range []float64{math.NaN(), math.Inf(1)} {
			x := testSignal(12)
			x[bad] = v
			ys := map[string][]float64{}
			var err error
			ys["Filter"], err = Filter(c.b, []float64{1}, x)
			require.NoError(t, err)
			ys["SOSFilter"], err = SOSFilter(SOS{{c.b[0], c.b[1], c.b[2], 1, 0, 0}}, x)
			require.NoError(t, err)
			clean, err := Filter(c.b, []float64{1}, testSignal(12))
			require.NoError(t, err)
			for name, y := range ys {
				for i := range y {
					if slices.Contains(c.affected, i) {
						assert.False(t, isFinite(y[i]), "%s %s, %g, output %d", name, c.name, v, i)
						continue
					}
					assert.Equal(t, clean[i], y[i], "%s %s, %g, output %d", name, c.name, v, i)
				}
			}
		}
	}
}

// isFinite reports whether v is neither NaN nor infinite
func isFinite(v float64) bool {
	return !math.IsNaN(v) && !math.IsInf(v, 0)
}
