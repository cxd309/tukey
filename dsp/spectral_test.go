package dsp

import (
	"math"
	"math/cmplx"
	"math/rand/v2"
	"slices"
	"testing"

	"github.com/cxd309/tukey/fft"
	"github.com/cxd309/tukey/window"
	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

// ramp is 1, 2, ... n, whose segments' sums are easy to check by hand.
func ramp(n int) (x []float64) {
	x = make([]float64, n)
	for i := range x {
		x[i] = float64(i + 1)
	}
	return
}

// hann is a periodic Hann window of length n, the usual window for spectral analysis.
func hann(t *testing.T, n int) (w []float64) {
	t.Helper()
	w, err := window.Hann(n, window.Periodic())
	require.NoError(t, err)
	return
}

// noise is n samples of standard normal noise, the same every run.
func noise(n int, seed uint64) (x []float64) {
	rng := rand.New(rand.NewPCG(seed, seed+1))
	x = make([]float64, n)
	for i := range x {
		x[i] = rng.NormFloat64()
	}
	return
}

// assertComplexNear asserts got is within tol of want, element by element,
// relative to want's largest magnitude.
func assertComplexNear(t *testing.T, want, got []complex128, tol float64) {
	t.Helper()
	if !assert.Len(t, got, len(want)) {
		return
	}
	var scale float64
	for _, w := range want {
		scale = max(scale, cmplx.Abs(w))
	}
	for i := range want {
		assert.LessOrEqual(t, cmplx.Abs(got[i]-want[i]), tol*scale, "[%d]: got %v, want %v", i, got[i], want[i])
	}
}

// Segments are cut by hand here: a 4-sample boxcar over 1 to 10 at 2 Hz.
// The default overlap is half the window, so segments start every 2 samples,
// and each one's DC bin is its sum.
func TestSTFTSegments(t *testing.T) {
	boxcar := []float64{1, 1, 1, 1}
	cases := []struct {
		name  string
		x     []float64
		opts  []SpectralOption
		times []float64
		sums  []float64
	}{
		{"default overlap", ramp(10), nil, []float64{1, 2, 3, 4}, []float64{10, 18, 26, 34}},
		{"partial last segment dropped", ramp(11), nil, []float64{1, 2, 3, 4}, []float64{10, 18, 26, 34}},
		{"overlap 1", ramp(10), []SpectralOption{Overlap(1)}, []float64{1, 2.5, 4}, []float64{10, 22, 34}},
		{"no overlap", ramp(10), []SpectralOption{Overlap(0)}, []float64{1, 3}, []float64{10, 26}},
		{"one segment", ramp(4), nil, []float64{1}, []float64{10}},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			tf, err := STFT(c.x, 2, boxcar, c.opts...)
			require.NoError(t, err)
			assert.Equal(t, c.times, tf.Times)
			assert.Equal(t, []float64{0, 0.5, 1}, tf.Frequencies) // NFFT is the window's length by default
			require.Len(t, tf.Values, len(c.sums))
			for i, sum := range c.sums {
				require.Len(t, tf.Values[i], 3)
				assert.Equal(t, complex(sum, 0), tf.Values[i][0], "segment %d's DC bin", i)
			}
		})
	}
}

// An odd window's default overlap rounds down, so segments start every 3 samples.
func TestSTFTDefaultOverlapOddWindow(t *testing.T) {
	tf, err := STFT(ramp(11), 1, []float64{1, 1, 1, 1, 1})
	require.NoError(t, err)
	assert.Equal(t, []float64{2.5, 5.5, 8.5}, tf.Times)
}

// At the FFT's own bins, AtFrequencies' direct sums agree with the FFT.
func TestSTFTAtFrequenciesMatchesBins(t *testing.T) {
	x, win := noise(600, 7), hann(t, 64)
	const fs = 100.0
	bins, err := STFT(x, fs, win)
	require.NoError(t, err)
	at, err := STFT(x, fs, win, AtFrequencies(bins.Frequencies))
	require.NoError(t, err)
	assert.Equal(t, bins.Frequencies, at.Frequencies)
	assert.Equal(t, bins.Times, at.Times)
	require.Len(t, at.Values, len(bins.Values))
	for i := range bins.Values {
		assertComplexNear(t, bins.Values[i], at.Values[i], 1e-13)
	}
}

// A real signal's negative frequencies are the conjugates of its positive ones,
// and two-sided frequencies are ordered as fft.Freq orders them.
func TestSTFTTwoSided(t *testing.T) {
	x, win := noise(300, 8), hann(t, 32)
	const fs = 50.0
	cases := []struct {
		name string
		nfft int
	}{
		{"even", 32},
		{"odd", 33},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			two, err := STFT(x, fs, win, NFFT(c.nfft), TwoSided())
			require.NoError(t, err)
			one, err := STFT(x, fs, win, NFFT(c.nfft))
			require.NoError(t, err)
			assert.Equal(t, fft.Freq(c.nfft, 1/fs), two.Frequencies)
			for i, spectrum := range two.Values {
				require.Len(t, spectrum, c.nfft)
				mirrored := make([]complex128, c.nfft)
				for k := range spectrum {
					mirrored[k] = cmplx.Conj(spectrum[(c.nfft-k)%c.nfft])
				}
				assertComplexNear(t, spectrum, mirrored, 1e-13)
				assertComplexNear(t, one.Values[i], spectrum[:len(one.Values[i])], 1e-13)
			}
		})
	}
}

// With each segment's mean removed, an offset added to the whole signal makes
// no difference. Without, it lands in the DC bin, scaled by the window's sum.
func TestSTFTSegmentDetrend(t *testing.T) {
	x, win := noise(256, 9), hann(t, 32)
	offset := make([]float64, len(x))
	for i, v := range x {
		offset[i] = v + 5
	}
	plain, err := STFT(x, 1, win, SegmentDetrend(DetrendConstant))
	require.NoError(t, err)
	shifted, err := STFT(offset, 1, win, SegmentDetrend(DetrendConstant))
	require.NoError(t, err)
	for i := range plain.Values {
		assertComplexNear(t, plain.Values[i], shifted.Values[i], 1e-13)
	}

	plain, err = STFT(x, 1, win)
	require.NoError(t, err)
	shifted, err = STFT(offset, 1, win)
	require.NoError(t, err)
	var windowSum float64
	for _, w := range win {
		windowSum += w
	}
	assert.InDelta(t, 5*windowSum, real(shifted.Values[0][0]-plain.Values[0][0]), 1e-12)
}

// STFT changes neither its inputs nor the frequencies given to AtFrequencies,
// even if the caller changes that slice after making the option.
func TestSTFTDoesNotShareInputs(t *testing.T) {
	x, win := noise(128, 10), hann(t, 16)
	frequencies := []float64{1, 2, 3}
	xBefore, winBefore := slices.Clone(x), slices.Clone(win)

	option := AtFrequencies(frequencies)
	frequencies[0] = 99
	tf, err := STFT(x, 10, win, option)
	require.NoError(t, err)
	assert.Equal(t, []float64{1, 2, 3}, tf.Frequencies)
	assert.Equal(t, xBefore, x)
	assert.Equal(t, winBefore, win)

	tf.Frequencies[1] = 99 // a second call isn't affected by changes to the first's result
	again, err := STFT(x, 10, win, option)
	require.NoError(t, err)
	assert.Equal(t, []float64{1, 2, 3}, again.Frequencies)
}

func TestSTFTInvalid(t *testing.T) {
	x, win := ramp(16), []float64{1, 1, 1, 1}
	cases := []struct {
		name   string
		x      []float64
		fs     float64
		window []float64
		opts   []SpectralOption
		want   error
	}{
		{"empty window", x, 1, nil, nil, ErrInvalidLength},
		{"window longer than x", x[:3], 1, win, nil, ErrInvalidLength},
		{"NFFT shorter than window", x, 1, win, []SpectralOption{NFFT(3)}, ErrInvalidLength},
		{"no frequencies", x, 1, win, []SpectralOption{AtFrequencies(nil)}, ErrInvalidLength},
		{"zero fs", x, 0, win, nil, ErrInvalidParameter},
		{"negative fs", x, -1, win, nil, ErrInvalidParameter},
		{"NaN fs", x, math.NaN(), win, nil, ErrInvalidParameter},
		{"negative overlap", x, 1, win, []SpectralOption{Overlap(-1)}, ErrInvalidParameter},
		{"overlap of the whole window", x, 1, win, []SpectralOption{Overlap(4)}, ErrInvalidParameter},
		{"AtFrequencies with NFFT", x, 1, win, []SpectralOption{AtFrequencies([]float64{0.1}), NFFT(8)}, ErrInvalidParameter},
		{"AtFrequencies two-sided", x, 1, win, []SpectralOption{AtFrequencies([]float64{0.1}), TwoSided()}, ErrInvalidParameter},
		{"unknown detrend mode", x, 1, win, []SpectralOption{SegmentDetrend(DetrendMode(99))}, ErrInvalidMode},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			_, err := STFT(c.x, c.fs, c.window, c.opts...)
			assert.ErrorIs(t, err, c.want)
		})
	}
}
