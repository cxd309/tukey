package floats

import "slices"

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
