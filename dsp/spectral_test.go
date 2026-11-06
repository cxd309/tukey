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

// Parseval: density summed over frequency is the windowed signal's power,
// sum((w·x)²) / sum(w²), so a boxcar's is x's mean square. One-sided, at even
// and odd lengths, as doubling every bin but DC and Nyquist keeps the total.
func TestPeriodogramParseval(t *testing.T) {
	const fs = 50.0
	for _, n := range []int{256, 255} {
		x := noise(n, 11)
		boxcar := make([]float64, n)
		for i := range boxcar {
			boxcar[i] = 1
		}
		for name, win := range map[string][]float64{"boxcar": boxcar, "hann": hann(t, n)} {
			t.Run(name+"/"+map[bool]string{true: "even", false: "odd"}[n%2 == 0], func(t *testing.T) {
				spectrum, err := Periodogram(x, fs, win)
				require.NoError(t, err)
				var total, power, squares float64
				for _, p := range spectrum.Power {
					total += p * fs / float64(n)
				}
				for i, v := range x {
					power += (win[i] * v) * (win[i] * v)
					squares += win[i] * win[i]
				}
				assert.InDelta(t, power/squares, total, 1e-12*power/squares)
			})
		}
	}
}

// With ScalingSpectrum, a sine on a bin reads its power, A²/2, at its frequency.
func TestPeriodogramSineSpectrum(t *testing.T) {
	const n, fs, amplitude = 64, 64.0, 3.0
	x := make([]float64, n)
	boxcar := make([]float64, n)
	for i := range x {
		x[i] = amplitude * math.Sin(2*math.Pi*8*float64(i)/fs) // 8 Hz is bin 8
		boxcar[i] = 1
	}
	spectrum, err := Periodogram(x, fs, boxcar, Scaling(ScalingSpectrum))
	require.NoError(t, err)
	assert.Equal(t, 8.0, spectrum.Frequencies[8])
	assert.InDelta(t, amplitude*amplitude/2, spectrum.Power[8], 1e-12)
}

// One-sided power is two-sided power folded over: each positive frequency plus
// its negative mirror, with DC, and Nyquist for an even NFFT, left as they are.
func TestSpectrogramOneSidedFoldsTwoSided(t *testing.T) {
	x, win := noise(400, 12), hann(t, 32)
	for _, nfft := range []int{32, 33} {
		t.Run(map[bool]string{true: "even", false: "odd"}[nfft%2 == 0], func(t *testing.T) {
			one, err := Spectrogram(x, 1, win, NFFT(nfft))
			require.NoError(t, err)
			two, err := Spectrogram(x, 1, win, NFFT(nfft), TwoSided())
			require.NoError(t, err)
			require.Len(t, one.Values, len(two.Values))
			for i, folded := range two.Values {
				want := make([]float64, nfft/2+1)
				want[0] = folded[0]
				for k := 1; k < len(want); k++ {
					want[k] = folded[k] + folded[nfft-k]
					if nfft%2 == 0 && k == nfft/2 {
						want[k] = folded[k] // Nyquist is its own mirror
					}
				}
				assert.InDeltaSlice(t, want, one.Values[i], 1e-12*slices.Max(want), "segment %d", i)
			}
		})
	}
}

// Over a single segment there's nothing to average, so Welch is that segment's
// periodogram, with either average: the median of one value has no bias.
func TestWelchOneSegmentIsPeriodogram(t *testing.T) {
	x, win := noise(64, 13), hann(t, 64)
	want, err := Periodogram(x, 10, win)
	require.NoError(t, err)
	for _, average := range []SpectralAverage{AverageMean, AverageMedian} {
		got, err := Welch(x, 10, win, Average(average))
		require.NoError(t, err)
		assert.Equal(t, want, got)
	}
}

// The median of n exponentially distributed values falls below their mean by
// 1 − 1/2 + 1/3 − … to the 1/n term for odd n, approaching ln 2 as n grows.
func TestMedianBias(t *testing.T) {
	cases := []struct {
		n    int
		want float64
	}{
		{1, 1},
		{2, 1},
		{3, 1 - 1.0/2 + 1.0/3},
		{4, 1 - 1.0/2 + 1.0/3},
		{5, 1 - 1.0/2 + 1.0/3 - 1.0/4 + 1.0/5},
	}
	for _, c := range cases {
		assert.InDelta(t, c.want, medianBias(c.n), 1e-15, "n=%d", c.n)
	}
	assert.InDelta(t, math.Ln2, medianBias(10001), 1e-4)
}

func TestPowerSpectraInvalid(t *testing.T) {
	x, win := ramp(16), []float64{1, 1, 1, 1}
	cases := []struct {
		name string
		call func() error
		want error
	}{
		{"Periodogram window shorter than x", func() error { _, err := Periodogram(x, 1, win); return err }, ErrInvalidLength},
		{"Overlap on Periodogram", func() error { _, err := Periodogram(x[:4], 1, win, Overlap(1)); return err }, ErrInvalidParameter},
		{"Average on Periodogram", func() error { _, err := Periodogram(x[:4], 1, win, Average(AverageMedian)); return err }, ErrInvalidParameter},
		{"Average on Spectrogram", func() error { _, err := Spectrogram(x, 1, win, Average(AverageMedian)); return err }, ErrInvalidParameter},
		{"AtFrequencies on Spectrogram", func() error { _, err := Spectrogram(x, 1, win, AtFrequencies([]float64{0.1})); return err }, ErrInvalidParameter},
		{"AtFrequencies on Welch", func() error { _, err := Welch(x, 1, win, AtFrequencies([]float64{0.1})); return err }, ErrInvalidParameter},
		{"Scaling on STFT", func() error { _, err := STFT(x, 1, win, Scaling(ScalingSpectrum)); return err }, ErrInvalidParameter},
		{"Average on STFT", func() error { _, err := STFT(x, 1, win, Average(AverageMedian)); return err }, ErrInvalidParameter},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			assert.ErrorIs(t, c.call(), c.want)
		})
	}
}
