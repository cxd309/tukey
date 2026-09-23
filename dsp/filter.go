package dsp

import "fmt"

// digitalFilter is a validated IIR filter in transfter-function form
// normalised so a[0] == 1, with b and a zero-padded to the same length
type digitalFilter struct {
	b, a []float64
}

// newDigitalFilter validates caller-supplied coefficients and normalises them
// copies are taken, so the caller's slices are never modified
func newDigitalFilter(b, a []float64) (f digitalFilter, err error) {
	if len(b) == 0 || len(a) == 0 {
		return f, fmt.Errorf("dsp: b and a must be non-empty, got len(b)=%d len(a)=%d", len(b), len(a))
	}
	if a[0] == 0 {
		return f, fmt.Errorf("dsp: a[0] must be non-zero")
	}
	n := max(len(b), len(a))
	f.b = make([]float64, n)
	f.a = make([]float64, n)
	for i, v := range b {
		f.b[i] = v / a[0]
	}
	for i, v := range a {
		f.a[i] = v / a[0]
	}
	return
}

// apply runs the filter over x in transposed direct form II, starting from
// the internal state initial (length n-1), or at rest when initial is nil
// equivalent to scipy.signal.lfilter
func (f digitalFilter) apply(x, initial []float64) (y []float64) {
	n := len(f.a)
	// one spare slot on the end stays 0, so the last state update needs no special case
	state := make([]float64, n)
	copy(state, initial)

	y = make([]float64, len(x))
	for i, xi := range x {
		yi := f.b[0]*xi + state[0]
		for k := range n - 1 {
			state[k] = f.b[k+1]*xi - f.a[k+1]*yi + state[k+1]
		}
		y[i] = yi
	}
	return
}

func Filter(b, a, x []float64) (y []float64, err error) {
	f, err := newDigitalFilter(b, a)
	if err != nil {
		return nil, err
	}
	y = f.apply(x, nil)
	return
}
