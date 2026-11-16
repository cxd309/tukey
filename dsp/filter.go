package dsp

import (
	"fmt"
	"math"

	"github.com/cxd309/tukey/internal/floats"
)

// Filter applies the IIR filter, b, a to x causally (one direction only)
// start from rest, so the output is phase-shifted; use FiltFilt for zero phase
//
// a must describe a stable filter (all poles strictly inside the unit circle);
// this isn't checked, matching scipy.signal.lfilter and MATLAB's filter,
// and an unstable filter produces output that grows without bound
//
// returns ErrInvalidCoefficients if b or a is empty, a[0] is zero, or any coefficient is NaN or infinite
//
// equivalent to scipy.signal.lfilter(b, a, x) and MATLAB's filter(b, a, x)
func Filter(b, a, x []float64) (y []float64, err error) {
	f, err := newDigitalFilter(b, a)
	if err != nil {
		return nil, err
	}
	y, _ = f.apply(x, nil)
	return
}

// FilterState is Filter starting from state zi, of length max(len(a), len(b)) - 1,
// and also returns the final state zf. Passing each block's zf as the next block's zi
// filters a long signal in blocks exactly as filtering it in one go would
//
// returns ErrInvalidCoefficients as Filter does, and ErrInvalidState if zi has the wrong length or isn't finite
//
// equivalent to scipy.signal.lfilter(b, a, x, zi=zi)
func FilterState(b, a, x, zi []float64) (y, zf []float64, err error) {
	f, err := newDigitalFilter(b, a)
	if err != nil {
		return nil, nil, err
	}
	if len(zi) != len(f.a)-1 {
		return nil, nil, fmt.Errorf("%w: zi has %d values; this filter needs max(len(a), len(b)) - 1 = %d",
			ErrInvalidState, len(zi), len(f.a)-1)
	}
	if !allFinite(zi) {
		return nil, nil, fmt.Errorf("%w: zi must be finite", ErrInvalidState)
	}
	y, zf = f.apply(x, zi)
	return
}

// FilterZi returns the state the filter settles into after a long unit step input
// scaled by a signal's first sample, it starts FilterState without a startup transient
//
// returns ErrInvalidCoefficients as Filter does,
// ErrNoSteadyState if the filter has a pole at DC (sum(a) == 0)
//
// equivalent to scipy.signal.lfilter_zi(b, a)
func FilterZi(b, a []float64) (zi []float64, err error) {
	f, err := newDigitalFilter(b, a)
	if err != nil {
		return nil, err
	}
	zi, err = f.steadyState()
	return
}

// SOSFilter applies the filter sos to x causally, section by section, starting from rest
// stays accurate at high orders and narrow bands, where Filter's (b, a) form loses precision
//
// returns ErrInvalidSOS if sos has no sections, any section's a0 isn't 1, or any coefficient is NaN or infinite
//
// equivalent to scipy.signal.sosfilt(sos, x) and MATLAB's sosfilt(sos, x)
func SOSFilter(sos SOS, x []float64) (y []float64, err error) {
	c, err := sos.toCascade()
	if err != nil {
		return nil, err
	}
	y, _ = c.apply(x, nil)
	return
}

// SOSFilterState is SOSFilter starting from state zi, one pair per section, and also returns the final state zf
// passing each block's zf as the next block's zi filters a long signal in blocks exactly as filtering it in one go would
//
// returns ErrInvalidSOS as SOSFilter does, and ErrInvalidState if zi doesn't have one pair per section or isn't finite
//
// equivalent to scipy.signal.sosfilt(sos, x, zi=zi)
func SOSFilterState(sos SOS, x []float64, zi [][2]float64) (y []float64, zf [][2]float64, err error) {
	c, err := sos.toCascade()
	if err != nil {
		return nil, nil, err
	}
	if len(zi) != len(sos) {
		return nil, nil, fmt.Errorf("%w: zi has %d pairs but sos has %d sections", ErrInvalidState, len(zi), len(sos))
	}
	for i, pair := range zi {
		if !allFinite(pair[:]) {
			return nil, nil, fmt.Errorf("%w: zi[%d] must be finite", ErrInvalidState, i)
		}
	}
	y, final := c.apply(x, toCascadeState(zi))
	zf = fromCascadeState(final)
	return
}

// SOSFilterZi returns the state each section settles into after a long unit step input
// scaled by a signal's first sample, it starts SOSFilterState without a startup transient
//
// returns ErrInvalidSOS as SOSFilter does, and ErrNoSteadyState if a section has a pole at DC
//
// equivalent to scipy.signal.sosfilt_zi(sos)
func SOSFilterZi(sos SOS) (zi [][2]float64, err error) {
	c, err := sos.toCascade()
	if err != nil {
		return nil, err
	}
	state, err := c.steadyState()
	if err != nil {
		return nil, err
	}
	zi = fromCascadeState(state)
	return
}

// digitalFilter is a validated IIR filter in transfter-function form
// normalised so a[0] == 1, with b and a zero-padded to the same length
type digitalFilter struct {
	b, a []float64
}

// newDigitalFilter validates caller-supplied coefficients and normalises them
// copies are taken, so the caller's slices are never modified
func newDigitalFilter(b, a []float64) (f digitalFilter, err error) {
	if err = validateCoefficients(b, a); err != nil {
		return f, err
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

// validateCoefficients checks that b and a describe a filter: both non-empty, and a[0] non-zero
func validateCoefficients(b, a []float64) (err error) {
	if len(b) == 0 || len(a) == 0 {
		return fmt.Errorf("%w: b and a must be non-empty, got len(b)=%d len(a)=%d", ErrInvalidCoefficients, len(b), len(a))
	}
	if a[0] == 0 {
		return fmt.Errorf("%w: a[0] must be non-zero", ErrInvalidCoefficients)
	}
	if !allFinite(b) || !allFinite(a) {
		return fmt.Errorf("%w: b and a must be finite", ErrInvalidCoefficients)
	}
	return
}

// apply runs the filter over x in transposed direct form II, starting from
// the internal state initial (length n-1), or at rest when initial is nil
// returns the output and the state it finishes in
// equivalent to scipy.signal.lfilter
func (f digitalFilter) apply(x, initial []float64) (y, final []float64) {
	n := len(f.a)
	// one spare slot on the end stays 0, so the last state update needs no special case
	state := make([]float64, n)
	copy(state, initial)

	y = make([]float64, len(x))
	for i, xi := range x {
		yi := f.b[0]*xi + state[0]
		if math.IsNaN(yi) {
			yi = term(f.b[0], xi) + state[0]
		}
		for k := range n - 1 {
			next := f.b[k+1]*xi - f.a[k+1]*yi + state[k+1]
			if math.IsNaN(next) {
				next = term(f.b[k+1], xi) - term(f.a[k+1], yi) + state[k+1]
			}
			state[k] = next
		}
		y[i] = yi
	}
	final = state[:n-1]
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
	sumB, sumA := f.sums()
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

// sums returns sum(b) and sum(a): the numerator and denominator at z = 1, which is DC.
// They're summed with compensation, as a sharp filter's coefficients nearly cancel
func (f digitalFilter) sums() (sumB, sumA float64) {
	var b, a floats.Sum
	for i := range f.a {
		b.Add(f.b[i])
		a.Add(f.a[i])
	}
	sumB, sumA = b.Total(), a.Total()
	return
}

// cascade is a filter as sections applied in turn, each section's output feeding the next
// a single (b, a) filter is a cascade of one, and an SOS is a cascade of second-order sections
type cascade []digitalFilter

// apply runs x through each section in turn, starting each from its state in initial
// every section at rest when initial is nil
// returns the output and each section's final state
// equivalent to scipy.signal.sosfilter, or lfilter for a single section
func (c cascade) apply(x []float64, initial [][]float64) (y []float64, final [][]float64) {
	y = x
	final = make([][]float64, len(c))
	for i, f := range c {
		var zi []float64
		if initial != nil {
			zi = initial[i]
		}
		y, final[i] = f.apply(y, zi)
	}
	return
}

// steadyState returns each section's state after a long unit step into the cascade
// each section's own steady state, scaled by the DC gain of the sections before it
// equivalent to scipy.signal.sosfilt_zi, or lfilter_zi for a single section
func (c cascade) steadyState() (zi [][]float64, err error) {
	zi = make([][]float64, len(c))
	scale := 1.0
	for i, f := range c {
		sectionZi, err := f.steadyState()
		if err != nil {
			return nil, err
		}
		zi[i] = floats.Scaled(sectionZi, scale)
		sumB, sumA := f.sums()
		scale *= sumB / sumA
	}
	return
}

// scaledState returns a copy of every section's state multiplied by s
func scaledState(zi [][]float64, s float64) (scaled [][]float64) {
	scaled = make([][]float64, len(zi))
	for i := range zi {
		scaled[i] = floats.Scaled(zi[i], s)
	}
	return
}

// toCascadeState converts per-section states from SOS's fixed pairs to cascade's slices
func toCascadeState(zi [][2]float64) (state [][]float64) {
	state = make([][]float64, len(zi))
	for i := range zi {
		state[i] = []float64{zi[i][0], zi[i][1]}
	}
	return
}

// fromCascadeState converts per-section states from cascade's slices to SOS's fixed pairs
func fromCascadeState(state [][]float64) (zi [][2]float64) {
	zi = make([][2]float64, len(state))
	for i := range state {
		zi[i] = [2]float64{state[i][0], state[i][1]}
	}
	return
}

// term is c·v, or 0 when c is exactly 0, so a coefficient of 0 drops a NaN or infinite v,
// which 0·v would turn into NaN.
// apply uses it only to recompute a NaN, so every finite result is unchanged.
func term(c, v float64) (t float64) {
	if c == 0 {
		return 0
	}
	t = c * v
	return
}
