package waveform

import (
	"fmt"
	"math"

	"github.com/cxd309/tukey/internal/special"
)

// SweepPoly returns cos(2π·φ(t) + phi) at each sample time in t,
// a cosine whose frequency, the derivative of φ, is the polynomial poly at t,
// its coefficients highest power first.
//
//	x, err := waveform.SweepPoly(t, []float64{0.025, -0.36, 1.25, 2}) // 0.025t³ - 0.36t² + 1.25t + 2 Hz
//
// The phase φ is the polynomial's integral from 0, in cycles, computed in double-double arithmetic,
// then its whole cycles removed before the cosine, same as Chirp.
// A NaN or infinite time gives NaN at that sample.
//
// Units:
//   - t is in seconds, or any time unit
//   - poly is in hertz, its coefficient of t^k in hertz per second^k
//   - x is dimensionless, between -1 and 1
//
// Returns:
//   - ErrInvalidParameter for an empty poly, or a NaN or infinite coefficient or Phase.
//   - ErrInvalidOption for an option other than Phase.
//
// Equivalent to:
//   - scipy.signal.sweep_poly(t, poly, phi)
//   - MATLAB has no equivalent
func SweepPoly(t, poly []float64, opts ...Option) (x []float64, err error) {
	cfg, err := newConfig("SweepPoly", opts, optionSet{optionPhase: true})
	if err != nil {
		return nil, err
	}
	if err = checkSweepPoly(poly, cfg.phase); err != nil {
		return nil, err
	}
	integral := integrate(poly)
	x = cosines(t, cfg.phase, func(t float64) (cycles special.DoubleDouble) {
		// Horner's rule, then a last multiply by t, as the integral's constant term is 0
		at := special.Exact(t)
		cycles = integral[0]
		for _, c := range integral[1:] {
			cycles = cycles.Mul(at).Add(c)
		}
		cycles = cycles.Mul(at)
		return
	})
	return
}

// checkSweepPoly returns ErrInvalidParameter for an empty poly, or a non-finite coefficient or phi.
func checkSweepPoly(poly []float64, phi float64) (err error) {
	if len(poly) == 0 {
		return fmt.Errorf("%w: SweepPoly's poly is empty, needs at least 1 coefficient", ErrInvalidParameter)
	}
	for i, c := range poly {
		if math.IsNaN(c) || math.IsInf(c, 0) {
			return fmt.Errorf("%w: SweepPoly's poly[%d] is %g, must be finite", ErrInvalidParameter, i, c)
		}
	}
	if math.IsNaN(phi) || math.IsInf(phi, 0) {
		return fmt.Errorf("%w: SweepPoly's Phase is %g, must be finite", ErrInvalidParameter, phi)
	}
	return
}

// integrate is poly's integral, highest power first, without its constant term:
// the coefficient of t^k becomes that of t^(k+1), divided by k+1
func integrate(poly []float64) (integral []special.DoubleDouble) {
	integral = make([]special.DoubleDouble, len(poly))
	for i, c := range poly {
		power := len(poly) - i // of t, once integrated
		integral[i] = special.Exact(c).Div(special.Exact(float64(power)))
	}
	return
}
