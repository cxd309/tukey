package window

import (
	"fmt"
	"math"
)

// Gaussian returns the n-point Gaussian window with standard deviation std, in samples
//
// returns ErrInvalidLength if n < 0, and ErrInvalidParameter unless std > 0
//
// equivalent to scipy.signal.windows.gaussian(n, std),
// equivalent to MATLAB sgausswin(n, alpha) with alpha = (n-1)/(2*std)
func Gaussian(n int, std float64, opts ...Option) (w []float64, err error) {
	if err = positive("Gaussian", "std", std); err != nil {
		return nil, err
	}
	w, err = newWindow("Gaussian", n, opts, func(k, m int) float64 {
		x := float64(k) - float64(m-1)/2
		return math.Exp(-x * x / (2 * math.Pow(std, 2)))
	})

	return
}

// GeneralGaussian returns the n-point generalised Gaussian window,
// exp(-0.5 |x/sig|^(2p)) for x samples from the centre
// p = 1 is Gaussian with std sig, larger p flattens the top
//
// returns ErrInvalidLength if n < 0, and ErrInvalidParameter unless p and sig are > 0
//
// equivalent to scipy.signal.windows.general_gaussian(n, p, sig); MATLAB has none
func GeneralGaussian(n int, p, sig float64, opts ...Option) (w []float64, err error) {
	if err = positive("GeneralGaussian", "p", p); err != nil {
		return nil, err
	}
	if err = positive("GeneralGaussian", "sig", sig); err != nil {
		return nil, err
	}
	w, err = newWindow("GeneralGaussian", n, opts, func(k, m int) float64 {
		x := float64(k) - float64(m-1)/2
		return math.Exp(-0.5 * math.Pow(math.Abs(x/sig), 2*p))
	})
	return
}

// Exponential returns the n-point exponential (Poisson) window,
// exp(-|x|/tau) for x samples from its peak, the middle unless Center moves it.
// The usual force and response window in impact testing
//
// returns ErrInvalidLength if n < 0, and ErrInvalidParameter unless tau > 0,
// or if Center is given without Periodic or is NaN
//
// equivalent to scipy.signal.windows.exponential(n, center, tau); MATLAB has none
func Exponential(n int, tau float64, opts ...Option) (w []float64, err error) {
	cfg, err := newConfig("Exponential", opts, optionSet{optionPeriodic: true, optionCenter: true})
	if err != nil {
		return nil, err
	}
	if err = positive("Exponential", "tau", tau); err != nil {
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
		return nil, fmt.Errorf("%w: Exponential's Center needs Periodic: a symmetric window is centred", ErrInvalidParameter)
	case math.IsNaN(cfg.center):
		return nil, fmt.Errorf("%w: Exponential's Center is NaN", ErrInvalidParameter)
	}
	// moved off the middle the window isn't symmetric, so every sample is computed
	w, err = sampled(n, func(i int) float64 {
		return math.Exp(-math.Abs(float64(i)-cfg.center) / tau)
	})
	return
}
