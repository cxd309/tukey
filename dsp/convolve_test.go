package dsp

import (
	"math"
	"slices"
	"testing"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

// convolveMethods are every ConvolveMethod, named for subtests.
var convolveMethods = []struct {
	name   string
	method ConvolveMethod
}{
	{"auto", ConvolveAuto},
	{"direct", ConvolveDirect},
	{"FFT", ConvolveFFT},
}

// assertConvolved checks got against want: exactly for the direct method,
// whose small integer sums are exact, and to rounding for the others.
func assertConvolved(t *testing.T, method ConvolveMethod, want, got []float64, msg string) {
	t.Helper()
	if method == ConvolveDirect {
		assert.Equal(t, want, got, msg)
		return
	}
	assert.InDeltaSlice(t, want, got, 1e-12, msg)
}

func TestConvolve(t *testing.T) {
	cases := []struct {
		name              string
		a, v              []float64
		full, same, valid []float64
	}{
		// same starts at (3 − 1)/2 = 1
		{"odd kernel", []float64{1, 2, 3}, []float64{0, 1, 0.5}, []float64{0, 1, 2.5, 4, 1.5}, []float64{1, 2.5, 4}, []float64{2.5}},
		// same starts at (2 − 1)/2 = 0, where MATLAB's starts at 1
		{"even kernel", []float64{1, 2, 3, 4}, []float64{1, 1}, []float64{1, 3, 5, 7, 4}, []float64{1, 3, 5, 7}, []float64{3, 5, 7}},
		// valid swaps the inputs, where MATLAB's is empty
		{"kernel longer", []float64{1, 1}, []float64{1, 2, 3}, []float64{1, 3, 5, 3}, []float64{3, 5}, []float64{3, 5}},
		{"one sample each", []float64{3}, []float64{-2}, []float64{-6}, []float64{-6}, []float64{-6}},
	}
	for _, c := range cases {
		for _, m := range convolveMethods {
			t.Run(c.name+", "+m.name, func(t *testing.T) {
				by := ConvolveBy(m.method)
				full, err := Convolve(c.a, c.v, by)
				require.NoError(t, err)
				assertConvolved(t, m.method, c.full, full, "full")
				same, err := Convolve(c.a, c.v, ConvolveSame(), by)
				require.NoError(t, err)
				assertConvolved(t, m.method, c.same, same, "same")
				valid, err := Convolve(c.a, c.v, ConvolveValid(), by)
				require.NoError(t, err)
				assertConvolved(t, m.method, c.valid, valid, "valid")
			})
		}
	}
}

// Convolution is commutative, so full and valid don't depend on the order of
// the inputs. Same does, as it keeps the first input's length.
func TestConvolveIsCommutative(t *testing.T) {
	a, v := []float64{1, -2, 3, 0.5, 4}, []float64{2, 0.25, -1}
	for _, opts := range [][]ConvolveOption{nil, {ConvolveValid()}} {
		av, err := Convolve(a, v, append(opts, ConvolveBy(ConvolveDirect))...)
		require.NoError(t, err)
		va, err := Convolve(v, a, append(opts, ConvolveBy(ConvolveDirect))...)
		require.NoError(t, err)
		assert.Equal(t, av, va)
	}
}

// The direct method keeps each product exactly and sums with compensation, so
// a small output between large ones that cancel keeps every digit.
func TestConvolveDirectKeepsSmallOutputs(t *testing.T) {
	a := []float64{1e8, 1e-3, -1e8}
	full, err := Convolve(a, []float64{1, 1, 1}, ConvolveBy(ConvolveDirect))
	require.NoError(t, err)
	assert.Equal(t, 1e-3, full[2])
}

// An infinite or NaN sample affects only the outputs it's summed into, and an
// infinite product isn't turned into NaN by its rounding error.
func TestConvolveDirectNonFinite(t *testing.T) {
	inf := math.Inf(1)
	full, err := Convolve([]float64{1, inf, 2}, []float64{1, 1}, ConvolveBy(ConvolveDirect))
	require.NoError(t, err)
	assert.Equal(t, []float64{1, inf, inf, 2}, full)

	full, err = Convolve([]float64{1, nan, 2}, []float64{1, 1}, ConvolveBy(ConvolveDirect))
	require.NoError(t, err)
	assert.Equal(t, 1.0, full[0])
	assert.True(t, math.IsNaN(full[1]))
	assert.True(t, math.IsNaN(full[2]))
	assert.Equal(t, 2.0, full[3])
}

// Auto picks direct for short inputs and FFT for long ones, and returns
// exactly what the picked method does.
func TestConvolveAuto(t *testing.T) {
	assert.Equal(t, ConvolveDirect, fasterMethod(5, 4))
	assert.Equal(t, ConvolveFFT, fasterMethod(4000, 4000))
	for _, n := range []int{5, 4000} {
		a, v := ramp(n), ramp(n-1)
		auto, err := Convolve(a, v)
		require.NoError(t, err)
		picked, err := Convolve(a, v, ConvolveBy(fasterMethod(n, n-1)))
		require.NoError(t, err)
		assert.Equal(t, picked, auto, "%d samples", n)
	}
}

// Same and valid are parts of the full result, sliced so that appending to
// one reallocates rather than writing over the outputs after it.
func TestConvolvePartsHaveNoSpareCapacity(t *testing.T) {
	a, v := ramp(10), ramp(4)
	for _, opt := range []ConvolveOption{ConvolveSame(), ConvolveValid()} {
		c, err := Convolve(a, v, opt)
		require.NoError(t, err)
		assert.Equal(t, len(c), cap(c))
		corr, err := Correlate(a, v, opt)
		require.NoError(t, err)
		assert.Equal(t, len(corr.Values), cap(corr.Values))
	}
}

func TestConvolveAndCorrelateLeaveInputsUnchanged(t *testing.T) {
	a, v := []float64{1, 2, 3}, []float64{4, 5}
	for _, m := range convolveMethods {
		_, err := Convolve(a, v, ConvolveBy(m.method))
		require.NoError(t, err)
		_, err = Correlate(a, v, ConvolveBy(m.method))
		require.NoError(t, err)
		assert.Equal(t, []float64{1, 2, 3}, a, m.name)
		assert.Equal(t, []float64{4, 5}, v, m.name)
	}
}

func TestConvolveAndCorrelateErrors(t *testing.T) {
	for _, in := range [][2][]float64{{nil, {1}}, {{1}, nil}, {{}, {}}} {
		_, err := Convolve(in[0], in[1])
		assert.ErrorIs(t, err, ErrInvalidLength, "Convolve of %d and %d samples", len(in[0]), len(in[1]))
		_, err = Correlate(in[0], in[1])
		assert.ErrorIs(t, err, ErrInvalidLength, "Correlate of %d and %d samples", len(in[0]), len(in[1]))
	}
	for _, method := range []ConvolveMethod{ConvolveAuto - 1, ConvolveFFT + 1} {
		_, err := Convolve([]float64{1}, []float64{1}, ConvolveBy(method))
		assert.ErrorIs(t, err, ErrInvalidMode, "Convolve, method %d", method)
		_, err = Correlate([]float64{1}, []float64{1}, ConvolveBy(method))
		assert.ErrorIs(t, err, ErrInvalidMode, "Correlate, method %d", method)
	}
}

func TestCorrelate(t *testing.T) {
	cases := []struct {
		name   string
		a, v   []float64
		opts   []ConvolveOption
		values []float64
		lags   []int
	}{
		// Σ a[n+k]·v[n]: lag −2 is a[0]·v[2], lag 2 is a[2]·v[0]
		{"full", []float64{1, 2, 3}, []float64{0, 1, 0.5}, nil, []float64{0.5, 2, 3.5, 3, 0}, []int{-2, -1, 0, 1, 2}},
		{"same", []float64{1, 2, 3}, []float64{0, 1, 0.5}, []ConvolveOption{ConvolveSame()}, []float64{2, 3.5, 3}, []int{-1, 0, 1}},
		{"valid", []float64{1, 2, 3}, []float64{0, 1, 0.5}, []ConvolveOption{ConvolveValid()}, []float64{3.5}, []int{0}},
		// an odd len(a) and an even len(v), where scipy's correlation_lags
		// labels same one lag high: same starts at (2 − 1)/2 = 0, lag −1
		{"same, odd and even", []float64{1, 2, 3}, []float64{1, 1}, []ConvolveOption{ConvolveSame()}, []float64{1, 3, 5}, []int{-1, 0, 1}},
		{"valid, kernel longer", []float64{1, 2}, []float64{1, 0, -1}, []ConvolveOption{ConvolveValid()}, []float64{-2, 1}, []int{-1, 0}},
	}
	for _, c := range cases {
		for _, m := range convolveMethods {
			t.Run(c.name+", "+m.name, func(t *testing.T) {
				corr, err := Correlate(c.a, c.v, append(c.opts, ConvolveBy(m.method))...)
				require.NoError(t, err)
				assertConvolved(t, m.method, c.values, corr.Values, "values")
				assert.Equal(t, c.lags, corr.Lags)
			})
		}
	}
}

// A received copy of a sent signal, delayed by d samples, correlates most
// strongly with it at lag d, in every mode.
func TestCorrelateFindsDelay(t *testing.T) {
	sent := []float64{1, -2, 3, 0.5, -1}
	const delay = 7
	received := make([]float64, 20)
	copy(received[delay:], sent)
	for _, opts := range [][]ConvolveOption{nil, {ConvolveSame()}, {ConvolveValid()}} {
		for _, m := range convolveMethods {
			corr, err := Correlate(received, sent, append(opts, ConvolveBy(m.method))...)
			require.NoError(t, err)
			i := slices.Index(corr.Values, slices.Max(corr.Values))
			assert.Equal(t, delay, corr.Lags[i], m.name)
		}
	}
}
