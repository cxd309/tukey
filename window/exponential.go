package window

import (
	"fmt"
	"math"
)

// Gaussian returns the n-point Gaussian window with standard deviation std:
//
//	w(x) = e^(-x²/(2·std²)), for x samples from the centre
//
// Units:
//   - n is a number of samples
//   - std is in samples
//   - w is dimensionless
//
// Returns:
//   - ErrInvalidLength for a negative n.
//   - ErrInvalidParameter for a std that isn't positive and finite.
//   - ErrInvalidOption for an option other than Periodic.
//
// Equivalent to:
//   - scipy.signal.windows.gaussian(n, std)
//   - MATLAB gausswin(n, (n-1)/(2*std)), which has no periodic form
func Gaussian(n int, std float64, opts ...Option) (w []float64, err error) {
	if err = positiveFinite("Gaussian", "std", std); err != nil {
		return nil, err
	}
	w, err = newWindow("Gaussian", n, opts, func(k, m int) float64 {
		x := float64(k) - float64(m-1)/2
		return math.Exp(-x * x / (2 * math.Pow(std, 2)))
	})
	return
}

// GeneralGaussian returns the n-point generalised Gaussian window:
//
//	w(x) = e^(-|x/sig|^(2p)/2), for x samples from the centre
//
// p = 1 is Gaussian with std sig, and a larger p flattens the top.
//
// Units:
//   - n is a number of samples
//   - p is dimensionless
//   - sig is in samples
//   - w is dimensionless
//
// Returns:
//   - ErrInvalidLength for a negative n.
//   - ErrInvalidParameter for a p or sig that isn't positive and finite.
//   - ErrInvalidOption for an option other than Periodic.
//
// Equivalent to:
//   - scipy.signal.windows.general_gaussian(n, p, sig)
//   - MATLAB has none
func GeneralGaussian(n int, p, sig float64, opts ...Option) (w []float64, err error) {
	if err = positiveFinite("GeneralGaussian", "p", p); err != nil {
		return nil, err
	}
	if err = positiveFinite("GeneralGaussian", "sig", sig); err != nil {
		return nil, err
	}
	w, err = newWindow("GeneralGaussian", n, opts, func(k, m int) float64 {
		x := float64(k) - float64(m-1)/2
		return math.Exp(-0.5 * math.Pow(math.Abs(x/sig), 2*p))
	})
	return
}

// Exponential returns the n-point exponential, or Poisson, window:
//
//	w(x) = e^(-|x|/tau), for x samples from its peak
//
// The peak is in the middle, unless Center moves it.
// It's the usual force and response window in impact testing.
//
//	w, err := window.Exponential(1024, 200, window.Periodic(), window.Center(0)) // decaying from the first sample
//
// Units:
//   - n is a number of samples
//   - tau is in samples, the distance over which w falls by a factor of e
//   - w is dimensionless
//
// Returns:
//   - ErrInvalidLength for a negative n.
//   - ErrInvalidParameter for a tau that isn't positive and finite, or a NaN or infinite Center.
//   - ErrInvalidOption for an option other than Periodic and Center, or Center without Periodic.
//
// Equivalent to:
//   - scipy.signal.windows.exponential(n, center, tau)
//   - MATLAB has none
func Exponential(n int, tau float64, opts ...Option) (w []float64, err error) {
	cfg, err := newConfig("Exponential", opts, optionSet{optionPeriodic: true, optionCenter: true})
	if err != nil {
		return nil, err
	}
	if err = positiveFinite("Exponential", "tau", tau); err != nil {
		return nil, err
	}
	if !cfg.given[optionCenter] {
		w, err = shape(n, cfg.periodic, func(k, m int) float64 {
			return math.Exp(-math.Abs(float64(k)-float64(m-1)/2) / tau)
		})
		return
	}
	switch {
	case !cfg.periodic:
		return nil, fmt.Errorf("%w: Exponential's Center needs Periodic: a symmetric window is centred", ErrInvalidOption)
	case math.IsNaN(cfg.center) || math.IsInf(cfg.center, 0):
		return nil, fmt.Errorf("%w: Exponential's Center %g, must be finite", ErrInvalidParameter, cfg.center)
	}
	// moved off the middle the window isn't symmetric, so every sample is computed
	w, err = sampled(n, func(i int) float64 {
		return math.Exp(-math.Abs(float64(i)-cfg.center) / tau)
	})
	return
}
