package dsp

import (
	"fmt"

	"github.com/cxd309/tukey/internal/poly"
)

// GroupDelay returns the group delay of the filter b/a in samples,
// at n frequencies w evenly spaced from 0 to just below π, as Freqz's:
// how many samples each frequency's envelope is delayed by, −d(phase)/dω.
//
//	w, delay, err := dsp.GroupDelay(b, a, 512)
//
// It's computed from the coefficients exactly, to rounding,
// not by differentiating the phase.
// A symmetric or antisymmetric b or a, such as a linear-phase FIR or a Butterworth numerator,
// contributes exactly half its order.
// Where b or a is ill-conditioned, its poles or zeros crowded against the unit circle,
// the point is evaluated in big.Float.
// At a zero on the unit circle, where the phase jumps by π, the delay is its limit from either side.
//
// A high-order filter's (b, a) is ill-conditioned,
// SOSGroupDelay on its sections is faster and gives the same answer.
//
// Returns:
//   - ErrInvalidCoefficients for coefficients Freqz rejects
//   - ErrInvalidLength if n < 0
//
// Equivalent to:
//   - scipy.signal.group_delay((b, a), n), which is 0 at a zero on the unit
//     circle and loses digits near one: 3% off for butter(8, 0.01)
//   - matlab [delay, w] = grpdelay(b, a, n), which also gives the limit
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

// SOSGroupDelay is GroupDelay for a filter given as second-order sections:
// the sum of each section's delay, which is better conditioned than the product's.
//
// Returns:
//   - ErrInvalidSOS for sections SOSFreqz rejects
//   - ErrInvalidLength if n < 0
//
// Equivalent to:
//   - scipy.signal.group_delay of the sections' product, which SciPy has no
//     sections form of
//   - matlab [delay, w] = grpdelay(sos, n)
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
