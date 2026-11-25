package waveform

import (
	"fmt"
	"math"

	"github.com/cxd309/tukey/internal/special"
)

// Sawtooth returns a sawtooth wave of period 2π at each sample time in t:
// rising from -1 to 1 over the first duty of each period, then falling back to -1 over the rest.
// Its duty is 1 by default, a rising ramp, or set with Duty, where 0 is a falling ramp and 0.5 a triangle.
//
//	x, err := waveform.Sawtooth(t, waveform.Duty(0.5)) // a triangle wave
//
// For a frequency f in hertz, pass `2* math.Pi* f* t`.
// t is reduced modulo the float64 2π, same as SciPy, so a t of 2π·k built from math.Pi is the start of a cycle.
// The ramps are computed in double-double from the exact remainder, so each sample is within an ULP of exact.
// A NaN or infinite time gives NaN at that sample.
//
// Units:
//   - t is in radians
//   - x is dimensionless, from -1 to 1
//
// Returns:
//   - ErrInvalidParameter for a Duty outside 0 and 1, or NaN
//   - ErrInvalidOption for an option other than Duty
//
// Equivalent to:
//   - scipy.signal.sawtooth(t, width), which returns NaN at width 1 for a tiny negative t,
//     where t% 2π rounds up to 2π
//   - matlab sawtooth(t, xmax), whose rem puts a t within round-off of 2π·k at the start of a cycle
func Sawtooth(t []float64, opts ...Option) (x []float64, err error) {
	cfg, err := newConfig("Sawtooth", opts, optionSet{optionDuty: true})
	if err != nil {
		return nil, err
	}
	duty := 1.0
	if cfg.given[optionDuty] {
		duty = cfg.duty
	}
	if !(duty >= 0 && duty <= 1) {
		return nil, fmt.Errorf("%w: Sawtooth's Duty is %g, must be between 0 and 1", ErrInvalidParameter, duty)
	}
	ramps := newSawtoothRamps(duty)
	x = make([]float64, len(t))
	for i, ti := range t {
		x[i] = ramps.at(ti)
	}
	return
}

// sawtoothRamps is where a sawtooth's rise ends, and how long its fall lasts, in radians
type sawtoothRamps struct {
	edge, fall special.DoubleDouble
}

// newSawtoothRamps is the ramps of a sawtooth rising for duty of each period of 2π
func newSawtoothRamps(duty float64) (ramps sawtoothRamps) {
	ramps = sawtoothRamps{
		edge: special.ProductOf(duty, 2*math.Pi),
		fall: special.SumOf(1, -duty).Mul(special.Exact(2 * math.Pi)),
	}
	return
}

// at is the sawtooth at t
func (ramps sawtoothRamps) at(t float64) (x float64) {
	if math.IsNaN(t) || math.IsInf(t, 0) {
		return math.NaN()
	}
	// mod.Math is exact, and the 2π added to a negative remainder is kept exact too,
	// so a tiny negative t stays just below 2π rather than rounded up to it
	r := special.Exact(math.Mod(t, 2*math.Pi))
	if r.Hi < 0 {
		r = special.SumOf(r.Hi, 2*math.Pi)
	}
	two, one := special.Exact(2), special.Exact(1)
	past := r.Sub(ramps.edge)
	if past.Hi < 0 {
		// rising: 2r/edge - 1, where edge is 0 only at duty 0, which never rises
		x = r.Mul(two).Div(ramps.edge).Sub(one).Float64()
		return
	}
	// falling: 1- 2(r- edge)/ fall, where fall is 0 only at duty 1, which never falls
	x = one.Sub(past.Mul(two).Div(ramps.fall)).Float64()
	return
}
