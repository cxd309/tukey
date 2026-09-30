package dsp

import (
	"math"
	"slices"
	"testing"

	"github.com/cxd309/tukey/internal/reference"
	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

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

// as for FiltFilt: steady-state starts and odd padding leave a constant untouched
func TestSOSFiltFiltPreservesConstant(t *testing.T) {
	sos := testSOS(t, 3, Lowpass(0.3))
	x := slices.Repeat([]float64{3}, 50)

	y, err := SOSFiltFilt(sos, x)
	require.NoError(t, err)
	reference.AssertClose(t, "y", y, x, reference.Tolerance{Rel: 1e-14, Scaled: 1e-14})
}

func TestSOSNtaps(t *testing.T) {
	cases := []struct {
		name  string
		order int
		want  int
	}{
		{"even order: 2 sections", 4, 5},
		{"odd order: 2 sections, one first-order", 3, 4},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			assert.Equal(t, c.want, testSOS(t, c.order, Lowpass(0.3)).ntaps())
		})
	}
}

func TestSOSFilteringRejectsInvalidInput(t *testing.T) {
	valid := SOS{{1, 0, 0, 1, -0.5, 0}}
	x := testSignal(50)

	_, err := SOSFilter(SOS{}, x)
	assert.ErrorIs(t, err, ErrInvalidSOS, "no sections")

	_, err = SOSFilter(SOS{{1, 0, 0, 2, -0.5, 0}}, x)
	assert.ErrorIs(t, err, ErrInvalidSOS, "a0 not 1")

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
