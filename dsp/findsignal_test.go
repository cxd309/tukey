package dsp

import (
	"math"
	"math/rand/v2"
	"slices"
	"testing"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

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
	methods := []struct {
		name   string
		method SearchMethod
	}{
		{"direct", SearchDirect},
		{"hybrid", SearchHybrid},
		{"fast", SearchFast},
	}
	for _, c := range findSignalBenchmarks() {
		for _, m := range methods {
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
