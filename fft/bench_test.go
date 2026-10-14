package fft

import (
	"fmt"
	"math"
	"testing"
)

// benchmarkLengths covers each algorithm path, at a small and large size
var benchmarkLengths = []struct {
	path string
	n    int
}{
	{"radix2", 1024},
	{"radix2", 65536},
	{"smooth", 1000},
	{"smooth", 48000}, // 2^7x3x5^3: one second at 48 kHz
	{"prime", 1009},   // a single Bluestein butterfly
	{"prime", 65537},
	{"largePrimeFactor", 4036}, // 4x1009: Bluestein butterflies inside mixed radix
}

func BenchmarkForward(b *testing.B) {
	for _, c := range benchmarkLengths {
		b.Run(fmt.Sprintf("%s/n=%d", c.path, c.n), func(b *testing.B) {
			x := randomSignal(c.n)
			b.ReportAllocs()
			b.ResetTimer()
			for range b.N {
				if _, err := Forward(x); err != nil {
					b.Fatal(err)
				}
			}
			reportPerNLogN(b, c.n)
		})
	}
}

func BenchmarkForwardReal(b *testing.B) {
	for _, c := range benchmarkLengths {
		b.Run(fmt.Sprintf("%s/n=%d", c.path, c.n), func(b *testing.B) {
			x := randomReal(c.n)
			b.ReportAllocs()
			b.ResetTimer()
			for range b.N {
				if _, err := ForwardReal(x); err != nil {
					b.Fatal(err)
				}
			}
			reportPerNLogN(b, c.n)
		})
	}
}

// reportPerNLogN adds the time per operation divided by n·log2(n)
// for an efficient FFT that should be about constant across lengths,
// so highlight an algorithm that is slow for its size
func reportPerNLogN(b *testing.B, n int) {
	perOp := float64(b.Elapsed().Nanoseconds()) / float64(b.N)
	b.ReportMetric(perOp/(float64(n)*math.Log2(float64(n))), "ns/nlogn")
}
