package dsp

import (
	"fmt"

	"github.com/cxd309/tukey/internal/poly"
)

// GroupDelay returns the group delay of the filter b/a at n frequencies w,
// evenly spaced from 0 up to, but not including, Nyquist, same as Freqz's:
// how long each frequency's envelope is delayed by, -dφ/dω.
// It's computed from the coefficients exactly, to rounding, not by differentiating the phase.
// A symmetric or antisymmetric b or a, such as a linear-phase FIR or a Butterworth numerator,
// contributes exactly half its order.
// Where b or a is ill-conditioned, its roots crowded against the unit circle,
// the point is evaluated in big.Float.
// At a zero on the unit circle, where the phase jumps by π, the delay is its limit from either side.
// A high-order filter's (b, a) is ill-conditioned, and SOSGroupDelay on its sections is faster and gives the same answer.
//
//	w, delay, err := dsp.GroupDelay(b, a, 512)
//
// Units:
//   - b and a are dimensionless
//   - w is in rad/sample, where π is Nyquist, and RadPerSampleToHz or RadPerSampleToNyquist converts it
//   - delay is in samples
//
// Returns:
//   - ErrInvalidCoefficients for an empty b or a, an a[0] of 0, or a NaN or infinite coefficient.
//   - ErrInvalidLength for a negative n.
//
// Equivalent to:
//   - scipy.signal.group_delay((b, a), n), which is 0 at a zero on the unit circle,
//     and loses digits near one, 3% off for butter(8, 0.01)
//   - MATLAB [delay, w] = grpdelay(b, a, n), which also gives the limit
func GroupDelay(b, a []float64, n int) (w, delay []float64, err error) {
	if err = validateCoefficients(b, a); err != nil {
		return nil, nil, err
	}
	if n < 0 {
		return nil, nil, fmt.Errorf("%w: n must be >= 0, got %d", ErrInvalidLength, n)
	}
	w = freqzPoints(n)
	delay = make([]float64, n)
	g := poly.NewDelayGrid(n)
	for i := range delay {
		delay[i] = g.Delay(b, i) - g.Delay(a, i)
	}
	return
}

// SOSGroupDelay returns the group delay of the filter sos at n frequencies w,
// evenly spaced from 0 up to, but not including, Nyquist, same as SOSFreqz's.
// It's the sum of each section's delay, computed exactly as GroupDelay's,
// which is better conditioned than the delay of the sections' product.
//
//	w, delay, err := dsp.SOSGroupDelay(sos, 512)
//
// Units:
//   - sos is dimensionless
//   - w is in rad/sample, where π is Nyquist, and RadPerSampleToHz or RadPerSampleToNyquist converts it
//   - delay is in samples
//
// Returns:
//   - ErrInvalidSOS for no sections, a section whose a0 isn't 1, or a NaN or infinite coefficient.
//   - ErrInvalidLength for a negative n.
//
// Equivalent to:
//   - scipy.signal.group_delay of the sections' product, as SciPy has no sections form
//   - MATLAB [delay, w] = grpdelay(sos, n)
func SOSGroupDelay(sos SOS, n int) (w, delay []float64, err error) {
	if err = sos.validate(); err != nil {
		return nil, nil, err
	}
	if n < 0 {
		return nil, nil, fmt.Errorf("%w: n must be >= 0, got %d", ErrInvalidLength, n)
	}
	w = freqzPoints(n)
	delay = make([]float64, n)
	g := poly.NewDelayGrid(n)
	for i := range delay {
		for _, section := range sos {
			delay[i] += g.Delay(section[:3], i) - g.Delay(section[3:], i)
		}
	}
	return
}
