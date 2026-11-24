package waveform

import (
	"fmt"
	"math"

	"github.com/cxd309/tukey/internal/special"
)

// Chirp returns cos(2π·φ(t) + phi) at each sample time in t,
// a cosine whose frequency, the derivative of φ, sweeps from f0 at t = 0 to f1 at t = t1.
// The sweep is Linear by default, or Quadratic, Logarithmic or Hyperbolic,
// and carries on past t1 and before 0.
//
//	t := make([]float64, 1000)
//	for i := range t {
//		t[i] = float64(i) / 1000 // one second at 1 kHz
//	}
//	x, err := waveform.Chirp(t, 10, 1, 200) // 10 Hz to 200 Hz
//
// The phase φ is in cycles and computed in double-double arithmetic,
// then its whole cycles removed before the cosine,
// so a sample is within a few ULP of exact however many cycles have passed.
// A NaN or infinite time gives NaN at that sample.
//
// Units:
//   - t and t1 are in seconds, or any time unit
//   - f0 and f1 are in hertz, or cycles per t's unit
//   - x is dimensionless, between -1 and 1
//
// Returns:
//   - ErrInvalidParameter for a NaN or infinite f0, t1, f1 or Phase, a t1 at or below 0,
//     an f0 or f1 of 0 or of opposite signs with Logarithmic, or an f0 or f1 of 0 with Hyperbolic.
//   - ErrInvalidOption for VertexAtT1 without Quadratic,
//     or an option other than Linear, Quadratic, Logarithmic, Hyperbolic, VertexAtT1 or Phase.
//
// Equivalent to:
//   - scipy.signal.chirp(t, f0, t1, f1, method, phi, vertex_zero)
//   - MATLAB chirp(t, f0, t1, f1, method, phi), which has no hyperbolic sweep,
//     and no quadratic with its vertex at t1
func Chirp(t []float64, f0, t1, f1 float64, opts ...Option) (x []float64, err error) {
	cfg, err := newConfig("Chirp", opts, optionSet{
		optionLinear: true, optionQuadratic: true, optionLogarithmic: true, optionHyperbolic: true,
		optionVertexAtT1: true, optionPhase: true,
	})
	if err != nil {
		return nil, err
	}
	if err = checkChirp(cfg, f0, t1, f1); err != nil {
		return nil, err
	}
	x = cosines(t, cfg.phase, chirpCycles(cfg, f0, t1, f1))
	return
}

// checkChirp returns ErrInvalidParameter or ErrInvalidOption for parameters the sweep in cfg can't take.
func checkChirp(cfg config, f0, t1, f1 float64) (err error) {
	for _, p := range []struct {
		name  string
		value float64
	}{{"f0", f0}, {"t1", t1}, {"f1", f1}, {"Phase", cfg.phase}} {
		if math.IsNaN(p.value) || math.IsInf(p.value, 0) {
			return fmt.Errorf("%w: Chirp's %s is %g, must be finite", ErrInvalidParameter, p.name, p.value)
		}
	}
	switch {
	case t1 <= 0:
		err = fmt.Errorf("%w: Chirp's t1 is %g, must be positive", ErrInvalidParameter, t1)
	case cfg.vertexAtT1 && cfg.sweep != sweepQuadratic:
		err = fmt.Errorf("%w: VertexAtT1 needs Quadratic", ErrInvalidOption)
	case cfg.sweep == sweepLogarithmic && (f0 == 0 || f1 == 0 || (f0 < 0) != (f1 < 0)):
		err = fmt.Errorf("%w: a Logarithmic Chirp's f0 and f1 are %g and %g, must be non-zero and of the same sign",
			ErrInvalidParameter, f0, f1)
	case cfg.sweep == sweepHyperbolic && (f0 == 0 || f1 == 0):
		err = fmt.Errorf("%w: a Hyperbolic Chirp's f0 and f1 are %g and %g, must be non-zero", ErrInvalidParameter, f0, f1)
	}
	return
}

// chirpCycles is the sweep's phase φ(t) in cycles, with φ(0) = 0, as a function of t.
// Each coefficient is computed once, in double-double.
func chirpCycles(cfg config, f0, t1, f1 float64) (cycles func(t float64) special.DoubleDouble) {
	rise := special.SumOf(f1, -f0)
	switch {
	case f0 == f1:
		// Logarithmic and Hyperbolic divide by f1 - f0 or its logarithm, and every sweep is constant here
		cycles = func(t float64) special.DoubleDouble { return special.ProductOf(f0, t) }

	case cfg.sweep == sweepQuadratic && cfg.vertexAtT1:
		// f1·t + (f1 - f0)·((t1 - t)³ - t1³)/(3·t1²)
		scale := rise.Div(special.ProductOf(t1, t1).Mul(special.Exact(3)))
		t1Cubed := cube(special.Exact(t1))
		cycles = func(t float64) special.DoubleDouble {
			return special.ProductOf(f1, t).Add(scale.Mul(cube(special.SumOf(t1, -t)).Sub(t1Cubed)))
		}

	case cfg.sweep == sweepQuadratic:
		// f0·t + (f1 - f0)·t³/(3·t1²)
		scale := rise.Div(special.ProductOf(t1, t1).Mul(special.Exact(3)))
		cycles = func(t float64) special.DoubleDouble {
			return special.ProductOf(f0, t).Add(scale.Mul(cube(special.Exact(t))))
		}

	case cfg.sweep == sweepLogarithmic:
		// t1·f0·((f1/f0)^(t/t1) - 1)/ln(f1/f0)
		logRatio := special.Log(special.Exact(f1).Div(special.Exact(f0)))
		rate := logRatio.Div(special.Exact(t1))
		scale := special.ProductOf(t1, f0).Div(logRatio)
		cycles = func(t float64) special.DoubleDouble {
			return scale.Mul(special.ExpM1(rate.Mul(special.Exact(t))))
		}

	case cfg.sweep == sweepHyperbolic:
		// -s·f0·ln|1 - t/s|, where s = f1·t1/(f1 - f0) is where the frequency is infinite
		singular := special.ProductOf(f1, t1).Div(rise)
		scale := singular.Mul(special.Exact(f0)).Neg()
		cycles = func(t float64) special.DoubleDouble {
			distance := special.Exact(1).Sub(special.Exact(t).Div(singular))
			if distance.Hi < 0 {
				distance = distance.Neg()
			}
			return scale.Mul(special.Log(distance))
		}

	default:
		// f0·t + (f1 - f0)·t²/(2·t1)
		scale := rise.Div(special.Exact(2 * t1))
		cycles = func(t float64) special.DoubleDouble {
			return special.ProductOf(f0, t).Add(scale.Mul(special.ProductOf(t, t)))
		}
	}
	return
}

// cube is d³.
func cube(d special.DoubleDouble) (cubed special.DoubleDouble) {
	cubed = d.Mul(d).Mul(d)
	return
}
