package dsp

import (
	"fmt"
	"math"

	"github.com/cxd309/tukey/internal/floats"
	"github.com/cxd309/tukey/internal/poly"
)

// SavGolOption configures SavGolFilter
type SavGolOption func(*savGolConfig)

// SavGolDeriv makes SavGolFilter return the order-th derivative of each fitted polynomial,
// for samples dx apart, rather than its value:
// a smoothed derivative.
// Order 0, the value, by default.
func SavGolDeriv(order int, dx float64) SavGolOption {
	return func(c *savGolConfig) { c.deriv, c.dx = order, dx }
}

// SavGolFilter smooths x with a Savitzky–Golay filter:
// each sample is replaced by a polynomial of degree polyOrder,
// fitted by least squares to the n samples centred on it, evaluated there.
// It smooths noise while keeping the height and width of peaks, which a moving average flattens.
// Within n/2 of the ends, one polynomial is fitted to the first or last n samples and evaluated at each.
// A window holding a NaN gives NaN.
//
// The fit is in polynomials orthogonal over the window's samples, built by Arnoldi iteration,
// so it stays accurate for long windows and high degrees, up to n − 1,
// where the fit passes through every sample.
// Each output is taken relative to its own sample, so a constant is kept exactly,
// and a small signal on a large offset keeps its digits.
//
// Returns:
//   - ErrInvalidLength for an n that isn't odd and positive
//   - ErrInvalidParameter for a polyOrder outside 0 to n − 1, or a
//     SavGolDeriv with a negative order or a dx that isn't positive and finite
//   - ErrSignalTooShort if x is shorter than n
//
// Equivalent to:
//   - scipy.signal.savgol_filter(x, n, polyOrder, deriv=order, delta=dx)
//   - matlab sgolayfilt(x, polyOrder, n), without SavGolDeriv
func SavGolFilter(x []float64, n, polyOrder int, opts ...SavGolOption) (y []float64, err error) {
	cfg := savGolConfig{dx: 1}
	for _, opt := range opts {
		opt(&cfg)
	}
	switch {
	case n < 1 || n%2 == 0:
		return nil, fmt.Errorf("%w: SavGolFilter needs an odd, positive window length, got %d", ErrInvalidLength, n)
	case polyOrder < 0 || polyOrder >= n:
		return nil, fmt.Errorf("%w: SavGolFilter polyOrder %d, must be from 0 to n − 1 = %d", ErrInvalidParameter, polyOrder, n-1)
	case cfg.deriv < 0:
		return nil, fmt.Errorf("%w: SavGolDeriv order %d, must be at least 0", ErrInvalidParameter, cfg.deriv)
	case !(cfg.dx > 0) || math.IsInf(cfg.dx, 1):
		return nil, fmt.Errorf("%w: SavGolDeriv dx %g, must be positive and finite", ErrInvalidParameter, cfg.dx)
	case len(x) < n:
		return nil, fmt.Errorf("%w: SavGolFilter needs at least n = %d samples, got %d", ErrSignalTooShort, n, len(x))
	}

	fit := poly.NewFit(n, polyOrder)
	// d/dt, for samples dx apart, is d/du over fit.Spacing()·dx
	scale := math.Pow(float64(fit.Spacing())*cfg.dx, float64(cfg.deriv))
	half := n / 2
	y = make([]float64, len(x))
	centre := fit.Weights(half, cfg.deriv)
	for i := half; i < len(x)-half; i++ {
		y[i] = fitted(centre, x[i-half:i+half+1], half, cfg.deriv) / scale
	}
	first, last := x[:n], x[len(x)-n:]
	for at := range half {
		y[at] = fitted(fit.Weights(at, cfg.deriv), first, at, cfg.deriv) / scale
		right := half + 1 + at
		y[len(x)-half+at] = fitted(fit.Weights(right, cfg.deriv), last, right, cfg.deriv) / scale
	}
	return
}

// savGolConfig is SavGolFilter's settings, built from SavGolOptions
type savGolConfig struct {
	deriv int
	dx    float64
}

// fitted is the fit's value, or deriv-th derivative, at window[at], from its weights w:
// c + Σ w_k·(x_k − c), with c = window[at],
// as the weights sum to 1, or to 0 for a derivative, which drops the c.
// Each x_k − c is recovered exactly by TwoSum, and every term summed with compensation,
// so a constant is kept exactly,
// and the rounding scales with how far the window's samples are from c rather than with their size:
// a small signal on a large offset keeps its digits.
// A non-finite c is left out, as 0, and a non-finite term is summed as it is,
// as there's no rounding error to recover from it.
func fitted(w, window []float64, at, deriv int) float64 {
	c := window[at]
	if math.IsInf(c, 0) || math.IsNaN(c) {
		c = 0
	}
	var sum floats.Sum
	if deriv == 0 {
		sum.Add(c)
	}
	for k, x := range window {
		d := x - c
		product := float64(w[k] * d)
		sum.Add(product)
		if math.IsInf(d, 0) || math.IsInf(product, 0) {
			continue
		}
		// TwoSum: x − c = d + e exactly
		shift := d - x
		e := (x - (d - shift)) + (-c - shift)
		sum.Add(math.FMA(w[k], d, -product))
		sum.Add(float64(w[k] * e))
	}
	return sum.Total()
}
