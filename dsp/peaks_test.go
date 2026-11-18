package dsp

import (
	"math"
	"testing"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

func TestFindPeaksLocalMaxima(t *testing.T) {
	cases := []struct {
		name string
		x    []float64
		want []int
	}{
		{"single", []float64{0, 1, 0}, []int{1}},
		{"ends are never peaks", []float64{3, 1, 2, 1, 3}, []int{2}},
		{"odd plateau, middle", []float64{0, 2, 2, 2, 0}, []int{2}},
		{"even plateau, middle rounded down", []float64{0, 2, 2, 0}, []int{1}},
		{"plateau to the end", []float64{0, 1, 2, 2}, []int{}},
		{"rising shoulder", []float64{0, 1, 1, 2, 0}, []int{3}},
		{"constant", []float64{1, 1, 1, 1}, []int{}},
		{"two samples", []float64{0, 1}, []int{}},
		{"empty", []float64{}, []int{}},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			peaks, err := FindPeaks(c.x)
			require.NoError(t, err)
			assert.Equal(t, c.want, peakIndices(peaks))
		})
	}
}

// Each peak is reported with the signal's value there.
func TestFindPeaksHeights(t *testing.T) {
	peaks, err := FindPeaks([]float64{0, 1, 0, 2.5, 0})
	require.NoError(t, err)
	assert.Equal(t, []Peak{{Index: 1, Height: 1}, {Index: 3, Height: 2.5}}, peaks)
}

// Each limit keeps the peaks that equal it.
func TestFindPeaksLimits(t *testing.T) {
	// peaks at 1, 3 and 5, of heights, thresholds and prominences 1, 2 and 3
	x := []float64{0, 1, 0, 2, 0, 3, 0}
	cases := []struct {
		name string
		opt  PeakOption
		want []int
	}{
		{"MinPeakHeight", MinPeakHeight(2), []int{3, 5}},
		{"PeakThreshold", PeakThreshold(2), []int{3, 5}},
		{"MinPeakProminence", MinPeakProminence(2), []int{3, 5}},
		{"MinPeakDistance at the gap", MinPeakDistance(2), []int{1, 3, 5}},
		{"MinPeakDistance over the gap", MinPeakDistance(3), []int{1, 5}},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			peaks, err := FindPeaks(x, c.opt)
			require.NoError(t, err)
			assert.Equal(t, c.want, peakIndices(peaks))
		})
	}
}

func TestFindPeaksPlateauFirst(t *testing.T) {
	peaks, err := FindPeaks([]float64{0, 2, 2, 2, 0}, PlateauFirst())
	require.NoError(t, err)
	assert.Equal(t, []int{1}, peakIndices(peaks))
}

func TestFindPeaksTies(t *testing.T) {
	x := []float64{0, 1, 0, 1, 0} // equal peaks 2 apart
	cases := []struct {
		name string
		ties PeakOption
		want []int
	}{
		{"TieLater", TieLater(), []int{3}},
		{"TieEarlier", TieEarlier(), []int{1}},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			peaks, err := FindPeaks(x, MinPeakDistance(3), c.ties)
			require.NoError(t, err)
			assert.Equal(t, c.want, peakIndices(peaks))
		})
	}
}

// A tall narrow spike at 2 next to a wide hump at 6.
// Applying distance first, the spike removes the hump and then fails the width itself,
// losing both. testdata's distance_vs_width_both is the same signal.
func TestFindPeaksDistanceLast(t *testing.T) {
	x := []float64{0, 0, 5, 0, 1, 2, 3, 2, 1, 0, 0}
	cases := []struct {
		name string
		opts []PeakOption
		want []int
	}{
		{"distance first", nil, []int{}},
		{"DistanceLast", []PeakOption{DistanceLast()}, []int{6}},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			opts := append([]PeakOption{MinPeakDistance(5), MinPeakWidth(2)}, c.opts...)
			peaks, err := FindPeaks(x, opts...)
			require.NoError(t, err)
			assert.Equal(t, c.want, peakIndices(peaks))
		})
	}
}

// MATLABPeaks sets each rule, and a later option overrides it.
func TestFindPeaksMATLABPeaksOverridden(t *testing.T) {
	x := []float64{0, 2, 2, 2, 0}
	cases := []struct {
		name string
		opts []PeakOption
		want []int
	}{
		{"MATLABPeaks", []PeakOption{MATLABPeaks()}, []int{1}},
		{"then PlateauMiddle", []PeakOption{MATLABPeaks(), PlateauMiddle()}, []int{2}},
		{"after PlateauMiddle", []PeakOption{PlateauMiddle(), MATLABPeaks()}, []int{1}},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			peaks, err := FindPeaks(x, c.opts...)
			require.NoError(t, err)
			assert.Equal(t, c.want, peakIndices(peaks))
		})
	}
}

// MATLAB errors where SciPy finds no peaks, and so does tukey.
func TestFindPeaksShortAndDistant(t *testing.T) {
	cases := []struct {
		name string
		x    []float64
		want []int
	}{
		{"empty", []float64{}, []int{}},
		{"one sample", []float64{1}, []int{}},
		{"two samples", []float64{0, 1}, []int{}},
		{"distance past the end", []float64{0, 1, 0}, []int{1}},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			peaks, err := FindPeaks(c.x, MinPeakDistance(len(c.x)+1), MATLABPeaks())
			require.NoError(t, err)
			assert.Equal(t, c.want, peakIndices(peaks))
		})
	}
}

func TestFindPeaksInvalidDistance(t *testing.T) {
	for _, d := range []int{0, -1} {
		_, err := FindPeaks([]float64{0, 1, 0}, MinPeakDistance(d))
		assert.ErrorIs(t, err, ErrInvalidParameter)
	}
}

func TestFindPeaksRejectsNaNLimits(t *testing.T) {
	nan := math.NaN()
	for name, opt := range map[string]PeakOption{
		"MinPeakHeight":     MinPeakHeight(nan),
		"PeakThreshold":     PeakThreshold(nan),
		"MinPeakProminence": MinPeakProminence(nan),
		"MinPeakWidth":      MinPeakWidth(nan),
	} {
		_, err := FindPeaks([]float64{0, 1, 0}, opt)
		assert.ErrorIs(t, err, ErrInvalidParameter, name)
	}
}

// An infinite limit is allowed:
// -Inf is no limit, and +Inf passes no peak.
func TestFindPeaksInfiniteLimits(t *testing.T) {
	x := []float64{0, 1, 0, 2, 0}
	inf := math.Inf(1)
	for name, opts := range map[string][2]PeakOption{
		"MinPeakHeight":     {MinPeakHeight(-inf), MinPeakHeight(inf)},
		"PeakThreshold":     {PeakThreshold(-inf), PeakThreshold(inf)},
		"MinPeakProminence": {MinPeakProminence(-inf), MinPeakProminence(inf)},
		"MinPeakWidth":      {MinPeakWidth(-inf), MinPeakWidth(inf)},
	} {
		all, err := FindPeaks(x, opts[0])
		require.NoError(t, err)
		assert.Equal(t, []int{1, 3}, peakIndices(all), "%s -Inf", name)
		none, err := FindPeaks(x, opts[1])
		require.NoError(t, err)
		assert.Empty(t, none, "%s +Inf", name)
	}
}

// A sample next to a NaN isn't known to be a peak, so it isn't one.
// A NaN reached while looking for a peak's base leaves its prominence and width unknown, so NaN,
// and a prominence or width limit drops it.
// A NaN beyond a higher sample isn't reached, so changes nothing.
func TestPeaksAroundNaN(t *testing.T) {
	nan := math.NaN()
	x := []float64{0, 3, 1, 4, 2, nan, 1, 5, 0}
	peaks, err := FindPeaks(x)
	require.NoError(t, err)
	assert.Equal(t, []int{1, 3, 7}, peakIndices(peaks))
	beside, err := FindPeaks([]float64{0, 2, nan, 0})
	require.NoError(t, err)
	assert.Empty(t, beside)

	prominences, err := PeakProminences(x, []int{1, 3, 7})
	require.NoError(t, err)
	assert.Equal(t, PeakProminence{Prominence: 2, LeftBase: 0, RightBase: 2}, prominences[0], "4 is reached before the NaN")
	assert.True(t, math.IsNaN(prominences[1].Prominence), "peak 3 reaches the NaN on its right")
	assert.True(t, math.IsNaN(prominences[2].Prominence), "peak 7 reaches the NaN on its left")

	widths, err := PeakWidths(x, []int{3}, 0.5)
	require.NoError(t, err)
	assert.True(t, math.IsNaN(widths[0].Width))

	kept, err := FindPeaks(x, MinPeakProminence(0.5))
	require.NoError(t, err)
	assert.Equal(t, []int{1}, peakIndices(kept))
	kept, err = FindPeaks(x, MinPeakWidth(0.5))
	require.NoError(t, err)
	assert.Equal(t, []int{1}, peakIndices(kept))
}

// Beside a -Inf the width's crossing is the next sample, its limit, where interpolating towards -Inf gives NaN.
func TestPeakWidthsBesideMinusInfinity(t *testing.T) {
	x := []float64{0, math.Inf(-1), 1, 2, 1, 0}
	widths, err := PeakWidths(x, []int{3}, 0.5)
	require.NoError(t, err)
	assert.Equal(t, 2.0, widths[0].Left)
	assert.Equal(t, 4.0, widths[0].Right)
}

func TestPeakProminences(t *testing.T) {
	x := []float64{0, 1, 0, 2, 0, 3, 0}
	got, err := PeakProminences(x, []int{1, 3, 5})
	require.NoError(t, err)
	want := []PeakProminence{
		{Prominence: 1, LeftBase: 0, RightBase: 2},
		{Prominence: 2, LeftBase: 2, RightBase: 4},
		{Prominence: 3, LeftBase: 4, RightBase: 6},
	}
	assert.Equal(t, want, got)
}

func TestPeakWidths(t *testing.T) {
	x := []float64{0, 1, 2, 1, 0} // a triangle of prominence 2
	cases := []struct {
		name      string
		relHeight float64
		want      PeakWidth
	}{
		{"on samples", 0.5, PeakWidth{Width: 2, Height: 1, Left: 1, Right: 3}},
		{"at the bases", 1, PeakWidth{Width: 4, Height: 0, Left: 0, Right: 4}},
		{"interpolated", 0.25, PeakWidth{Width: 1, Height: 1.5, Left: 1.5, Right: 2.5}},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			got, err := PeakWidths(x, []int{2}, c.relHeight)
			require.NoError(t, err)
			assert.Equal(t, []PeakWidth{c.want}, got)
		})
	}
}

func TestPeakProminencesAndWidthsInvalid(t *testing.T) {
	x := []float64{0, 1, 0}
	for _, peak := range []int{-1, len(x)} {
		_, err := PeakProminences(x, []int{peak})
		assert.ErrorIs(t, err, ErrInvalidParameter)
		_, err = PeakWidths(x, []int{peak}, 0.5)
		assert.ErrorIs(t, err, ErrInvalidParameter)
	}
	for _, relHeight := range []float64{-0.5, math.NaN(), math.Inf(1)} {
		_, err := PeakWidths(x, []int{1}, relHeight)
		assert.ErrorIs(t, err, ErrInvalidParameter)
	}
}

// peakIndices is the index of each of peaks.
func peakIndices(peaks []Peak) (indices []int) {
	indices = make([]int, len(peaks))
	for k, p := range peaks {
		indices[k] = p.Index
	}
	return
}
