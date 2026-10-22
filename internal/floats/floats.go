package floats

import (
	"math"
	"slices"
)

// Sum is a running sum with Neumaier's compensation:
// it keeps the low bits each addition round off, and adds them back in Total.
// A sum of many small values with a few large ones stays accurate,
// where a plain sum loses the small values' low bits.
// The zero Sum is empty.
type Sum struct {
	sum  float64
	lost float64 // the low bits rounded off so far
}

// Add adds v to the sum.
func (s *Sum) Add(v float64) {
	total := s.sum + v
	if math.Abs(s.sum) >= math.Abs(v) {
		s.lost += (s.sum - total) + v
	} else {
		s.lost += (v - total) + s.sum
	}
	s.sum = total
}

// Total is the sum of everything added, with the lost bits restored
func (s Sum) Total() (total float64) {
	total = s.sum + s.lost
	return
}

// Scaled returns a new slice holding v multiplied by s
func Scaled(v []float64, s float64) (out []float64) {
	out = make([]float64, len(v))
	for i := range v {
		out[i] = v[i] * s
	}
	return
}

// Reversed returns a new slice holding v in reverse order
func Reversed(v []float64) (r []float64) {
	r = slices.Clone(v)
	slices.Reverse(r)
	return
}
