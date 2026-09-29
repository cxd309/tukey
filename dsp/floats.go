package dsp

import "slices"

// scaled returns a new slice holding v multiplied by s
func scaled(v []float64, s float64) (out []float64) {
	out = make([]float64, len(v))
	for i := range v {
		out[i] = v[i] * s
	}
	return
}

// reversed returns a new slice holding v in reverse order
func reversed(v []float64) (r []float64) {
	r = slices.Clone(v)
	slices.Reverse(r)
	return
}
