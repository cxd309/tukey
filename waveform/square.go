package waveform

import (
	"fmt"
	"math"
)

// Square returns a square wave of period 2π at each sample time in t:
// 1 for the first duty of each period, and -1 for the rest.
// Its duty is 0.5 by default, or set with Duty.
//
//	x, err := waveform.Square(t, waveform.Duty(0.25)) // high for a quarter of each period
//
// For a frequency f in hertz, pass `2* math.Pi* f* t`.
// t is reduced modulo the float64 2π, same as SciPy, which is exact,
// so a t of 2π·k built from math.Pi is the start of a cycle, and gives 1.
// A NaN or infinite time gives NaN at that sample.
//
// Units:
//   - t is in radians
//   - x is dimensionless, 1 or -1
//
// Returns:
//   - ErrInvalidParameter for a Duty outside 0 to 1, or NaN.
//   - ErrInvalidOption for an option other than Duty.
//
// Equivalent to:
//   - scipy.signal.square(t, duty), which returns NaN for a duty outside 0 to -1,
//     and -1 for a NaN time
//   - MATLAB square(t, 100*duty), which takes duty as a percentage,
//     and whose mod puts a t within round-off of 2π·k at the start of a cycle
func Square(t []float64, opts ...Option) (x []float64, err error) {
	cfg, err := newConfig("Square", opts, optionSet{optionDuty: true})
	if err != nil {
		return nil, err
	}
	duty := 0.5
	if cfg.given[optionDuty] {
		duty = cfg.duty
	}
	if !(duty >= 0 && duty <= 1) {
		return nil, fmt.Errorf("%w: Square's Duty is %g, must be between 0 and 1", ErrInvalidParameter, duty)
	}
	edge := duty * 2 * math.Pi
	x = make([]float64, len(t))
	for i, ti := range t {
		x[i] = squareAt(ti, edge)
	}
	return
}

// squareAt is the square wave at t, which is 1 below edge in each period of 2π and -1 from edge on
func squareAt(t, edge float64) (x float64) {
	if math.IsNaN(t) || math.IsInf(t, 0) {
		return math.NaN()
	}
	// math.Mod is exact, and takes the sign of t, so a negative remainder is moved into a period
	r := math.Mod(t, 2*math.Pi)
	if r < 0 {
		r += 2 * math.Pi
	}
	x = -1
	if r < edge {
		x = 1
	}
	return
}
