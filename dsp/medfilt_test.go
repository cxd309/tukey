package dsp

import (
	"math"
	"math/big"
	"slices"
	"testing"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

// The cases use x = 5, 1, 4, 2, 8, 3, 9, whose truncated windows near the ends
// hold an even count for n = 5. MATLAB's medfilt1 gives the same.
func TestMedFilt(t *testing.T) {
	x := []float64{5, 1, 4, 2, 8, 3, 9}
	cases := []struct {
		name string
		n    int
		edge MedFiltEdge
		want []float64
	}{
		{"truncate n=3", 3, MedFiltTruncate, []float64{3, 4, 2, 4, 3, 8, 6}},
		{"truncate n=5", 5, MedFiltTruncate, []float64{4, 3, 4, 3, 4, 5.5, 8}},
		{"zero pad n=3", 3, MedFiltZeroPad, []float64{1, 4, 2, 4, 3, 8, 3}},
		{"zero pad n=5", 5, MedFiltZeroPad, []float64{1, 2, 4, 3, 4, 3, 3}},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			y, err := MedFilt(x, c.n, MedFiltEdges(c.edge))
			require.NoError(t, err)
			assert.Equal(t, c.want, y)
		})
	}
}

// The zero MedFiltEdge, and so no option, is MedFiltTruncate.
func TestMedFiltDefaultIsTruncate(t *testing.T) {
	var edge MedFiltEdge
	assert.Equal(t, MedFiltTruncate, edge)

	x := noise(50, 1)
	byDefault, err := MedFilt(x, 7)
	require.NoError(t, err)
	truncated, err := MedFilt(x, 7, MedFiltEdges(MedFiltTruncate))
	require.NoError(t, err)
	assert.Equal(t, truncated, byDefault)
}

// Against the median of each window sorted afresh, with its even-count
// midpoint exact, for every n up to past the signal's length, on both edges:
// on noise, on noise rounded to a few values, so windows hold ties, and with NaN.
func TestMedFiltMatchesSortingEachWindow(t *testing.T) {
	ties := noise(40, 2)
	for i := range ties {
		ties[i] = math.Round(ties[i])
	}
	withNaN := noise(40, 3)
	withNaN[0], withNaN[17], withNaN[18] = math.NaN(), math.NaN(), math.NaN()
	signals := map[string][]float64{"noise": noise(40, 1), "ties": ties, "nan": withNaN}

	for name, x := range signals {
		for _, edge := range []MedFiltEdge{MedFiltTruncate, MedFiltZeroPad} {
			for n := 1; n <= 45; n += 2 {
				y, err := MedFilt(x, n, MedFiltEdges(edge))
				require.NoError(t, err)
				for i := range x {
					want := medianBySorting(x, i, n, edge)
					if math.IsNaN(want) {
						assert.True(t, math.IsNaN(y[i]), "%s, edge %d, n %d, y[%d] = %v, want NaN", name, edge, n, i, y[i])
						continue
					}
					assert.Equal(t, want, y[i], "%s, edge %d, n %d, y[%d]", name, edge, n, i)
				}
			}
		}
	}
}

// medianBySorting is the median of the window of n samples centred on x[i],
// sorted afresh, with an even count's midpoint exact
func medianBySorting(x []float64, i, n int, edge MedFiltEdge) float64 {
	var window []float64
	for j := i - n/2; j <= i+n/2; j++ {
		switch {
		case j >= 0 && j < len(x):
			window = append(window, x[j])
		case edge == MedFiltZeroPad:
			window = append(window, 0)
		}
	}
	if slices.ContainsFunc(window, math.IsNaN) {
		return math.NaN()
	}
	slices.Sort(window)
	middle := len(window) / 2
	if len(window)%2 == 1 {
		return window[middle]
	}
	return exactMidpoint(window[middle-1], window[middle])
}

// exactMidpoint is (a + b)/2 for finite a and b, summed exactly in big.Float,
// halved exactly, and rounded once to float64, to nearest even
func exactMidpoint(a, b float64) float64 {
	sum := new(big.Float).SetPrec(2200).SetFloat64(a)
	sum.Add(sum, new(big.Float).SetFloat64(b))
	half := new(big.Float).SetMantExp(sum, -1)
	m, _ := half.Float64()
	return m
}

// midpoint is (a + b)/2 rounded once, including where a + b overflows and
// where the half is subnormal, which a/2 + b/2 rounds twice.
func TestMidpoint(t *testing.T) {
	tiny := math.SmallestNonzeroFloat64
	cases := []struct {
		name string
		a, b float64
		want float64
	}{
		{"exact", 1, 2, 1.5},
		{"overflowing sum", 0x1p1023, 0x1.8p1023, 0x1.4p1023},
		{"overflowing negative sum", -math.MaxFloat64, -math.MaxFloat64, -math.MaxFloat64},
		{"opposite extremes", -math.MaxFloat64, math.MaxFloat64, 0},
		// 1.5 subnormal ULPs, a tie, to even: 2; a/2 + b/2 gives 0 + 1
		{"subnormal tie", tiny, 2 * tiny, 2 * tiny},
		// 0.5 ULPs, a tie, to even: 0
		{"smallest subnormal and 0", 0, tiny, 0},
		{"infinity", 1, math.Inf(1), math.Inf(1)},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			assert.Equal(t, c.want, midpoint(c.a, c.b))
		})
	}
	assert.True(t, math.IsNaN(midpoint(math.Inf(-1), math.Inf(1))), "midpoint of −Inf and +Inf")

	// against the exact midpoint, across magnitudes and signs
	values := []float64{0, tiny, 3 * tiny, 0x1p-1022, 0x1.8p-1022, 0.1, 1, 1 + 0x1p-52, 3, 1e300, math.MaxFloat64}
	for _, a := range values {
		for _, b := range values {
			for _, signs := range [][2]float64{{1, 1}, {1, -1}, {-1, 1}, {-1, -1}} {
				x, y := signs[0]*a, signs[1]*b
				assert.Equal(t, exactMidpoint(x, y), midpoint(x, y), "midpoint(%v, %v)", x, y)
			}
		}
	}
}

// The median filter's purpose: spikes go, and a step stays sharp.
func TestMedFiltRemovesSpikesAndKeepsSteps(t *testing.T) {
	x := []float64{20, 20, 20, 95, 20, 20, 30, 30, -40, 30, 30, 30}
	want := []float64{20, 20, 20, 20, 20, 20, 30, 30, 30, 30, 30, 30}
	y, err := MedFilt(x, 3)
	require.NoError(t, err)
	assert.Equal(t, want, y)
}

// A NaN gives NaN wherever a window holds it, and nowhere else.
func TestMedFiltNaN(t *testing.T) {
	x := []float64{1, 2, 3, 4, math.NaN(), 6, 7, 8, 9}
	for _, edge := range []MedFiltEdge{MedFiltTruncate, MedFiltZeroPad} {
		y, err := MedFilt(x, 5, MedFiltEdges(edge))
		require.NoError(t, err)
		for i, v := range y {
			assert.Equal(t, i >= 2 && i <= 6, math.IsNaN(v), "edge %d, y[%d] = %v", edge, i, v)
		}
	}
}

// Infinities sort like any value: a spike of +Inf is removed, and an
// even-count window holding both infinities has no midpoint.
func TestMedFiltInfinity(t *testing.T) {
	inf := math.Inf(1)
	y, err := MedFilt([]float64{1, 2, inf, 4, 5}, 3)
	require.NoError(t, err)
	assert.Equal(t, []float64{1.5, 2, 4, 5, 4.5}, y)

	y, err = MedFilt([]float64{-inf, inf, 1}, 3)
	require.NoError(t, err)
	assert.True(t, math.IsNaN(y[0]), "y[0] = %v, want NaN", y[0])
	assert.Equal(t, []float64{1, inf}, y[1:])
}

// n = 1 is the identity, into a new slice.
func TestMedFiltOneSampleWindow(t *testing.T) {
	x := noise(20, 4)
	for _, edge := range []MedFiltEdge{MedFiltTruncate, MedFiltZeroPad} {
		y, err := MedFilt(x, 1, MedFiltEdges(edge))
		require.NoError(t, err)
		assert.Equal(t, x, y)
		y[0]++
		assert.NotEqual(t, x[0], y[0], "y shares x's storage")
	}
}

// A window longer than x: truncated, every window is the whole of x;
// zero-padded, every window holds all of x and n − 3 zeros, which outnumber
// the samples from n = 7. The padding is counted, not stored, so even the
// largest n works.
func TestMedFiltWindowLongerThanX(t *testing.T) {
	x := []float64{3, 1, 2}
	cases := []struct {
		n       int
		zeroPad []float64
	}{
		{5, []float64{1, 1, 1}}, // 0, 0, 1, 2, 3
		{7, []float64{0, 0, 0}},
		{math.MaxInt, []float64{0, 0, 0}},
	}
	for _, c := range cases {
		y, err := MedFilt(x, c.n)
		require.NoError(t, err)
		assert.Equal(t, []float64{2, 2, 2}, y, "n %d", c.n)

		y, err = MedFilt(x, c.n, MedFiltEdges(MedFiltZeroPad))
		require.NoError(t, err)
		assert.Equal(t, c.zeroPad, y, "n %d", c.n)
	}
}

func TestMedFiltLeavesXUnchanged(t *testing.T) {
	x := []float64{5, 1, 4, math.NaN(), 8, 3, 9}
	original := slices.Clone(x)
	_, err := MedFilt(x, 5)
	require.NoError(t, err)
	assert.True(t, slices.EqualFunc(original, x, func(a, b float64) bool {
		return a == b || math.IsNaN(a) && math.IsNaN(b)
	}))
}

func TestMedFiltEmpty(t *testing.T) {
	for _, edge := range []MedFiltEdge{MedFiltTruncate, MedFiltZeroPad} {
		y, err := MedFilt(nil, 3, MedFiltEdges(edge))
		require.NoError(t, err)
		assert.Empty(t, y)
	}
}

func TestMedFiltRejectsInvalidInput(t *testing.T) {
	x := []float64{1, 2, 3}
	for _, n := range []int{0, -1, -3, 2, 4, math.MinInt} {
		_, err := MedFilt(x, n)
		assert.ErrorIs(t, err, ErrInvalidLength, "n %d", n)
	}
	for _, edge := range []MedFiltEdge{-1, 2} {
		_, err := MedFilt(x, 3, MedFiltEdges(edge))
		assert.ErrorIs(t, err, ErrInvalidMode, "edge %d", edge)
	}
}
