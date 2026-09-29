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
		return f, fmt.Errorf("%w: b and a must be non-empty, got len(b)=%d len(a)=%d", ErrInvalidCoefficients, len(b), len(a))
	}
	if a[0] == 0 {
		return f, fmt.Errorf("%w: a[0] must be non-zero", ErrInvalidCoefficients)
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

// steadyState returns the state apply settles into after a long unit step input,
// so starting from it (scaled by the first sample) avoids a startup transient
//
// once settled, the output is the DC gain g = sum(b)/sum(a); substituting x=1, y=g
// into apply's state updates lets them be solved from the last one back:
// zi[k] = b[k+1] - a[k+1]*g + zi[k+1]
// equivalent to scipy.signal.lfilter_zi, which solves the same equations as a linear system
func (f digitalFilter) steadyState() (zi []float64, err error) {
	sumB, sumA := 0.0, 0.0
	for i := range f.a {
		sumB += f.b[i]
		sumA += f.a[i]
	}
	if sumA == 0 {
		return nil, fmt.Errorf("%w: it has a pole at z=1 (sum(a) == 0)", ErrNoSteadyState)
	}
	gain := sumB / sumA

	n := len(f.a)
	zi = make([]float64, n-1)
	next := 0.0 // zi[k+1]; zero past the last state, like apply's spare slot
	for k := n - 2; k >= 0; k-- {
		zi[k] = f.b[k+1] - f.a[k+1]*gain + next
		next = zi[k]
	}
	return
}

// Filter applies the IIR filter, b, a to x causally (one direction only)
// start from rest, so the output is phase-shifted; use FiltFilt for zero phase
//
// a must describe a stable filter (all poles strictly inside the unit circle);
// this isn't checked, matching scipy.signal.lfilter and MATLAB's filter,
// and an unstable filter produces output that grows without bound
//
// returns ErrInvalidCoefficients if b or a is empty, or a[0] is zero
//
// equivalent to scipy.signal.lfilter(b, a, x) and MATLAB's filter(b, a, x)
func Filter(b, a, x []float64) (y []float64, err error) {
	f, err := newDigitalFilter(b, a)
	if err != nil {
		return nil, err
	}
	y = f.apply(x, nil)
	return
}
