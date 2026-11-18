package dsp

import (
	"fmt"
	"math"

	"github.com/cxd309/tukey/internal/floats"
)

// Filter applies the filter b/a to x causally, in one direction, starting from rest:
//
//	a[0]·y[n] = Σ b[k]·x[n-k] - Σ a[k]·y[n-k], the second sum from k = 1
//
// The output is delayed by the filter's phase, and FiltFilt gives zero phase instead.
// The filter must be stable, with every pole inside the unit circle.
// That isn't checked, same as SciPy and MATLAB, and an unstable filter's output grows without bound.
// As the order rises or the band narrows, b and a lose digits, and SOSFilter keeps them.
//
//	y, err := dsp.Filter([]float64{1, 1}, []float64{2}, x) // the mean of each sample and the one before
//
// Units:
//   - b and a are dimensionless
//   - y is in x's unit
//
// Returns: ErrInvalidCoefficients for an empty b or a, an a[0] of 0, or a NaN or infinite coefficient.
//
// Equivalent to:
//   - scipy.signal.lfilter(b, a, x)
//   - MATLAB filter(b, a, x)
func Filter(b, a, x []float64) (y []float64, err error) {
	f, err := newDigitalFilter(b, a)
	if err != nil {
		return nil, err
	}
	y, _ = f.apply(x, nil)
	return
}

// FilterState applies the filter b/a to x causally, starting from the state zi,
// and returns the final state zf too.
// Passing each block's zf as the next block's zi filters a long signal in blocks,
// exactly as filtering it in one go would.
// zi has max(len(a), len(b)) - 1 values, all 0 for rest, or FilterZi's scaled by the first sample to start without a transient.
//
//	zi := make([]float64, max(len(a), len(b))-1)
//	for _, block := range blocks {
//		y, zi, err = dsp.FilterState(b, a, block, zi)
//	}
//
// Units:
//   - b and a are dimensionless
//   - zi, zf and y are in x's unit
//
// Returns:
//   - ErrInvalidCoefficients for an empty b or a, an a[0] of 0, or a NaN or infinite coefficient.
//   - ErrInvalidState for a zi of the wrong length, or one that isn't finite.
//
// Equivalent to:
//   - scipy.signal.lfilter(b, a, x, zi=zi)
//   - MATLAB [y, zf] = filter(b, a, x, zi)
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

// FilterZi returns the state the filter b/a settles into after a long unit step.
// Scaled by a signal's first sample, it starts FilterState without a startup transient:
//
//	zi, err := dsp.FilterZi(b, a)
//	for i := range zi {
//		zi[i] *= x[0]
//	}
//	y, _, err := dsp.FilterState(b, a, x, zi)
//
// Units: b, a and zi are dimensionless, zi per unit of step
//
// Returns:
//   - ErrInvalidCoefficients for an empty b or a, an a[0] of 0, or a NaN or infinite coefficient.
//   - ErrNoSteadyState for a filter with a pole at DC, where Σ a = 0.
//
// Equivalent to:
//   - scipy.signal.lfilter_zi(b, a)
//   - MATLAB has none, as filtfilt computes it internally
func FilterZi(b, a []float64) (zi []float64, err error) {
	f, err := newDigitalFilter(b, a)
	if err != nil {
		return nil, err
	}
	zi, err = f.steadyState()
	return
}

// SOSFilter applies the filter sos to x causally, section by section, starting from rest.
// It stays accurate at high orders and narrow bands, where Filter's (b, a) form loses digits.
// The filter must be stable, which isn't checked, same as SciPy and MATLAB.
//
//	f, err := dsp.Butter(8, dsp.Lowpass(0.1))
//	sos, err := f.SOS()
//	y, err := dsp.SOSFilter(sos, x)
//
// Units:
//   - sos is dimensionless
//   - y is in x's unit
//
// Returns: ErrInvalidSOS for no sections, a section whose a0 isn't 1, or a NaN or infinite coefficient.
//
// Equivalent to:
//   - scipy.signal.sosfilt(sos, x)
//   - MATLAB sosfilt(sos, x)
func SOSFilter(sos SOS, x []float64) (y []float64, err error) {
	c, err := sos.toCascade()
	if err != nil {
		return nil, err
	}
	y, _ = c.apply(x, nil)
	return
}

// SOSFilterState applies the filter sos to x causally, section by section, starting from the state zi,
// one pair per section, and returns the final state zf too.
// Passing each block's zf as the next block's zi filters a long signal in blocks,
// exactly as filtering it in one go would.
//
//	zi := make([][2]float64, len(sos))
//	for _, block := range blocks {
//		y, zi, err = dsp.SOSFilterState(sos, block, zi)
//	}
//
// Units:
//   - sos is dimensionless
//   - zi, zf and y are in x's unit
//
// Returns:
//   - ErrInvalidSOS for no sections, a section whose a0 isn't 1, or a NaN or infinite coefficient.
//   - ErrInvalidState for a zi without one pair per section, or one that isn't finite.
//
// Equivalent to:
//   - scipy.signal.sosfilt(sos, x, zi=zi)
//   - MATLAB has none, as its sosfilt takes no state
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

// SOSFilterZi returns the state each section of sos settles into after a long unit step.
// Scaled by a signal's first sample, it starts SOSFilterState without a startup transient.
//
// Units: sos and zi are dimensionless, zi per unit of step
//
// Returns:
//   - ErrInvalidSOS for no sections, a section whose a0 isn't 1, or a NaN or infinite coefficient.
//   - ErrNoSteadyState for a section with a pole at DC.
//
// Equivalent to:
//   - scipy.signal.sosfilt_zi(sos)
//   - MATLAB has none
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

// digitalFilter is a validated filter in transfer-function form,
// normalised so a[0] = 1, with b and a zero-padded to the same length.
type digitalFilter struct {
	b, a []float64
}

// newDigitalFilter validates a caller's coefficients and normalises them.
// It copies them, so the caller's slices are never modified.
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

// validateCoefficients returns ErrInvalidCoefficients unless b and a describe a filter:
// both non-empty, a[0] non-zero, and every coefficient finite.
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

// apply runs the filter over x in transposed direct form II, same as scipy.signal.lfilter,
// starting from the state initial, of length n-1, or at rest when initial is nil.
// It returns the output and the state it finishes in.
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

// steadyState returns the state apply settles into after a long unit step,
// so starting from it, scaled by the first sample, avoids a startup transient.
// Once settled, the output is the DC gain g = Σ b/Σ a.
// Substituting x = 1 and y = g into apply's state updates lets them be solved from the last one back:
//
//	zi[k] = b[k+1] - a[k+1]·g + zi[k+1]
//
// scipy.signal.lfilter_zi solves the same equations as a linear system.
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

// sums returns Σ b and Σ a, the numerator and denominator at z = 1, which is DC.
// They're summed with compensation, as a sharp filter's coefficients nearly cancel.
func (f digitalFilter) sums() (sumB, sumA float64) {
	var b, a floats.Sum
	for i := range f.a {
		b.Add(f.b[i])
		a.Add(f.a[i])
	}
	sumB, sumA = b.Total(), a.Total()
	return
}

// cascade is a filter as sections applied in turn, each section's output feeding the next.
// A single (b, a) filter is a cascade of one, and an SOS is a cascade of second-order sections.
type cascade []digitalFilter

// apply runs x through each section in turn, same as scipy.signal.sosfilt, or lfilter for one section,
// starting each from its state in initial, or at rest when initial is nil.
// It returns the output and each section's final state.
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

// steadyState returns each section's state after a long unit step into the cascade, same as scipy.signal.sosfilt_zi:
// each section's own steady state, scaled by the DC gain of the sections before it.
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

// scaledState returns a copy of every section's state, multiplied by s.
func scaledState(zi [][]float64, s float64) (scaled [][]float64) {
	scaled = make([][]float64, len(zi))
	for i := range zi {
		scaled[i] = floats.Scaled(zi[i], s)
	}
	return
}

// toCascadeState converts per-section states from SOS's fixed pairs to cascade's slices.
func toCascadeState(zi [][2]float64) (state [][]float64) {
	state = make([][]float64, len(zi))
	for i := range zi {
		state[i] = []float64{zi[i][0], zi[i][1]}
	}
	return
}

// fromCascadeState converts per-section states from cascade's slices to SOS's fixed pairs.
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
