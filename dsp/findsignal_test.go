package dsp

import (
	"math"
	"math/rand/v2"
	"slices"
	"testing"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

// searchMethods is every SearchMethod, named for subtests.
var searchMethods = []struct {
	name   string
	method SearchMethod
}{
	{"direct", SearchDirect},
	{"hybrid", SearchHybrid},
	{"fast", SearchFast},
}

func TestFindSignal(t *testing.T) {
	oneSample := []float64{3, 1, 4, 1, 5}
	target := 4.2
	diff := oneSample[2] - target // at run time: as a constant it would be exact, not float64's
	cases := []struct {
		name         string
		data, signal []float64
		want         SignalMatch
	}{
		{"exact match", []float64{0, 1, 5, 3, 1, 0}, []float64{5, 3, 1}, SignalMatch{Start: 2, Stop: 4, Distance: 0}},
		{"whole", []float64{2, 7, 1}, []float64{2, 7, 1}, SignalMatch{Start: 0, Stop: 2, Distance: 0}},
		{"at the end", []float64{0, 0, 0, 1, 2}, []float64{1, 2}, SignalMatch{Start: 3, Stop: 4, Distance: 0}},
		{"no exact match", []float64{0, 10, 0, 1, 3, 0}, []float64{1, 2}, SignalMatch{Start: 3, Stop: 4, Distance: 1}},
		{"one sample", oneSample, []float64{target}, SignalMatch{Start: 2, Stop: 2, Distance: diff * diff}},
		// every period is an exact match, and the earliest wins
		{"ties", []float64{0, 1, 2, 1, 0, 1, 2, 1, 0, 1, 2, 1}, []float64{0, 1, 2, 1}, SignalMatch{Start: 0, Stop: 3, Distance: 0}},
	}
	for _, c := range cases {
		for _, m := range searchMethods {
			t.Run(c.name+"/"+m.name, func(t *testing.T) {
				match, err := FindSignal(c.data, c.signal, SignalSearch(m.method))
				require.NoError(t, err)
				assert.Equal(t, c.want.Start, match.Start, "start")
				assert.Equal(t, c.want.Stop, match.Stop, "stop")
				if m.method == SearchFast {
					// rounds at the scale of the sums of squares, not of the distance
					assert.InDelta(t, c.want.Distance, match.Distance, 1e-13, "distance")
					return
				}
				assert.Equal(t, c.want.Distance, match.Distance, "distance")
			})
		}
	}
}

// The zero SearchMethod, and so FindSignal with no options, is SearchHybrid.
func TestFindSignalDefaultIsHybrid(t *testing.T) {
	data := []float64{3, 1, 4, 1, 5, 9, 2, 6, 5, 3}
	signal := []float64{5, 9, 2}
	want, err := FindSignal(data, signal, SignalSearch(SearchHybrid))
	require.NoError(t, err)
	got, err := FindSignal(data, signal)
	require.NoError(t, err)
	assert.Equal(t, want, got)
}

// Hybrid gives Direct's answer exactly, start and distance, including on
// periodic data, where many windows are nearly as close as the best. Fast
// finds the same start except where its rounding can't tell windows apart,
// and never gives a negative distance.
func TestFindSignalMethodsAgree(t *testing.T) {
	rng := rand.New(rand.NewPCG(5, 6))
	walk := make([]float64, 3000)
	for i := 1; i < len(walk); i++ {
		walk[i] = walk[i-1] + rng.NormFloat64()
	}
	noisyCopy := make([]float64, 300)
	for i := range noisyCopy {
		noisyCopy[i] = walk[1000+i] + 0.1*rng.NormFloat64()
	}
	periodic := make([]float64, 1000)
	for i := range periodic {
		periodic[i] = math.Sin(2 * math.Pi * float64(i) / 50)
	}
	cases := []struct {
		name         string
		data, signal []float64
		nearTies     bool // whether Fast's rounding can confuse the closest windows
	}{
		{"self match", walk, walk[1700:2000], false},
		{"noisy copy", walk, noisyCopy, false},
		{"short self match", walk, walk[400:500], false},
		{"periodic", periodic, periodic[:50], true},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			direct, err := FindSignal(c.data, c.signal, SignalSearch(SearchDirect))
			require.NoError(t, err)
			hybrid, err := FindSignal(c.data, c.signal, SignalSearch(SearchHybrid))
			require.NoError(t, err)
			fast, err := FindSignal(c.data, c.signal, SignalSearch(SearchFast))
			require.NoError(t, err)
			assert.Equal(t, direct, hybrid, "hybrid")
			assert.GreaterOrEqual(t, fast.Distance, 0.0, "fast distance")
			if !c.nearTies {
				assert.Equal(t, direct.Start, fast.Start, "fast start")
			}
		})
	}
}

// NaN in an early window doesn't stop Direct or Hybrid finding a later match.
// Fast's FFT spreads a NaN into every window, so it has no answer to check.
func TestFindSignalNaN(t *testing.T) {
	data := []float64{math.NaN(), 1, 2, 3, 0, 5, 6}
	signal := []float64{2, 3, 0}
	for _, m := range searchMethods[:2] {
		t.Run(m.name, func(t *testing.T) {
			match, err := FindSignal(data, signal, SignalSearch(m.method))
			require.NoError(t, err)
			assert.Equal(t, SignalMatch{Start: 2, Stop: 4, Distance: 0}, match)
		})
	}
}

func TestFindSignalInvalid(t *testing.T) {
	cases := []struct {
		name         string
		data, signal []float64
		opts         []FindSignalOption
		want         error
	}{
		{"empty signal", []float64{1, 2}, nil, nil, ErrInvalidLength},
		{"empty both", nil, nil, nil, ErrInvalidLength},
		{"signal longer than data", []float64{1, 2}, []float64{1, 2, 3}, nil, ErrInvalidLength},
		{"unknown method", []float64{1, 2}, []float64{1}, []FindSignalOption{SignalSearch(SearchMethod(99))}, ErrInvalidMode},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			_, err := FindSignal(c.data, c.signal, c.opts...)
			assert.ErrorIs(t, err, c.want)
		})
	}
}

// distanceProfile's error in every window is within a quarter of
// profileErrorBound, so the bound has 4x to spare.
func TestProfileErrorBound(t *testing.T) {
	rng := rand.New(rand.NewPCG(1, 2))
	walk := make([]float64, 5000)
	for i := 1; i < len(walk); i++ {
		walk[i] = walk[i-1] + rng.NormFloat64()
	}
	noise := make([]float64, 3000)
	for i := range noise {
		noise[i] = 1000 + rng.NormFloat64() // large values, small differences
	}
	periodic := make([]float64, 2000)
	for i := range periodic {
		periodic[i] = math.Sin(2 * math.Pi * float64(i) / 50)
	}
	cases := []struct {
		name         string
		data, signal []float64
	}{
		{"walk", walk, walk[2000:2500]},
		{"offset noise", noise, noise[100:400]},
		{"periodic", periodic, periodic[:50]},
		{"one sample", noise, noise[7:8]},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			profile, err := distanceProfile(c.data, c.signal)
			require.NoError(t, err)
			bound := profileErrorBound(c.data, c.signal)
			m := len(c.signal)
			var worst float64
			for i, p := range profile {
				worst = max(worst, math.Abs(p-windowDistance(c.data[i:i+m], c.signal)))
			}
			t.Logf("worst error %.2e, bound %.2e, %.0fx to spare", worst, bound, bound/worst)
			assert.LessOrEqual(t, 4*worst, bound)
		})
	}
}

// findSignalBenchmark is one FindSignal input to benchmark.
type findSignalBenchmark struct {
	name         string
	data, signal []float64
}

// findSignalBenchmarks are the app's scale, a short template in the same
// data, and periodic data, where many windows are nearly as close as the best.
func findSignalBenchmarks() (cases []findSignalBenchmark) {
	rng := rand.New(rand.NewPCG(3, 4))
	walk := make([]float64, 400000)
	for i := 1; i < len(walk); i++ {
		walk[i] = walk[i-1] + rng.NormFloat64()
	}
	periodic := make([]float64, 100000)
	for i := range periodic {
		periodic[i] = math.Sin(2 * math.Pi * float64(i) / 50)
	}
	cases = []findSignalBenchmark{
		{"app/n=400000/m=30000", walk, walk[200000:230000]},
		{"short/n=400000/m=100", walk, walk[123456:123556]},
		{"periodic/n=100000/m=500", periodic, periodic[:500]},
	}
	return
}

func BenchmarkFindSignal(b *testing.B) {
	for _, c := range findSignalBenchmarks() {
		for _, m := range searchMethods {
			b.Run(c.name+"/"+m.name, func(b *testing.B) {
				b.ReportAllocs()
				for range b.N {
					if _, err := FindSignal(c.data, c.signal, SignalSearch(m.method)); err != nil {
						b.Fatal(err)
					}
				}
				if m.method == SearchHybrid {
					b.StopTimer() // counting reruns the profile, which mustn't count towards ns/op
					b.ReportMetric(float64(hybridCandidates(b, c.data, c.signal)), "candidates")
				}
			})
		}
	}
}

// hybridCandidates is how many windows findSignalHybrid rescores exactly.
func hybridCandidates(b *testing.B, data, signal []float64) (count int) {
	profile, err := distanceProfile(data, signal)
	if err != nil {
		b.Fatal(err)
	}
	cutoff := slices.Min(profile) + 2*profileErrorBound(data, signal)
	for _, p := range profile {
		if p <= cutoff {
			count++
		}
	}
	return
}
