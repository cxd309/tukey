package dsp

import (
	"fmt"
	"math"

	"github.com/cxd309/tukey/internal/floats"
)

// SavGolOption configures SavGolFilter
type SavGolOption func(*savGolConfig)

// savGolConfig is SavGolFilter's settings, built from SavGolOptions
type savGolConfig struct {
	deriv int
	dx    float64
}

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

	fit := newPolyFit(n, polyOrder)
	// d/dt, for samples dx apart, is d/du over fit.spacing·dx
	scale := math.Pow(float64(fit.spacing)*cfg.dx, float64(cfg.deriv))
	half := n / 2
	y = make([]float64, len(x))
	centre := fit.weights(half, cfg.deriv)
	for i := half; i < len(x)-half; i++ {
		y[i] = fitted(centre, x[i-half:i+half+1], half, cfg.deriv) / scale
	}
	first, last := x[:n], x[len(x)-n:]
	for at := range half {
		y[at] = fitted(fit.weights(at, cfg.deriv), first, at, cfg.deriv) / scale
		right := half + 1 + at
		y[len(x)-half+at] = fitted(fit.weights(right, cfg.deriv), last, right, cfg.deriv) / scale
	}
	return
}

// polyFit is the polynomials q_0 to q_polyOrder orthogonal over a window's n positions,
// the samples' offsets from its centre divided by spacing, so from −1 to 1.
// The least-squares fit of degree polyOrder is then Σ q_j·⟨q_j, x⟩/‖q_j‖²,
// a projection onto each with no system to solve.
//
// They're built by Arnoldi iteration, as Vandermonde with Arnoldi
// (Brubeck, Nakatsukasa and Trefethen, 2021):
// u·q_j orthogonalised against every q_i before it, twice,
// which stays orthogonal at any degree,
// where the Gram polynomials' three-term recurrence fails beyond about 2√n. q_0 is 1,
// with ‖q_0‖² = n, so a moving average's weight is 1/n rounded once,
// rather than (1/√n)² rounded three times, every later q_j is normalised, ‖q_j‖² = 1.
type polyFit struct {
	n, polyOrder int
	spacing      int         // half the window, or 1 for n = 1
	q            [][]float64 // q[j][k] = q_j at position k
	norm         []float64   // ‖q_j‖², summed over the positions: n, then 1s
	h            [][]float64 // the recurrence: u·q_j = Σ_i≤j+1 h[i][j]·q_i
}

// newPolyFit is the orthogonal polynomials up to polyOrder over n positions
func newPolyFit(n, polyOrder int) (f polyFit) {
	f = polyFit{n: n, polyOrder: polyOrder, spacing: max(n/2, 1)}
	f.q = make([][]float64, polyOrder+1)
	f.norm = make([]float64, polyOrder+1)
	f.h = make([][]float64, polyOrder+1)
	for j := range f.h {
		f.h[j] = make([]float64, polyOrder)
	}
	f.q[0] = make([]float64, n)
	for k := range f.q[0] {
		f.q[0][k] = 1
	}
	f.norm[0] = float64(n)
	for j := range polyOrder {
		v := make([]float64, n)
		for k := range v {
			v[k] = f.position(k) * f.q[j][k]
		}
		for range 2 {
			for i := 0; i <= j; i++ {
				c := compensatedDot(f.q[i], v) / f.norm[i]
				f.h[i][j] += c
				for k := range v {
					v[k] -= float64(c * f.q[i][k])
				}
			}
		}
		length := math.Sqrt(compensatedDot(v, v))
		for k := range v {
			v[k] /= length
		}
		f.h[j+1][j] = length
		f.q[j+1] = v
		f.norm[j+1] = 1
	}
	return
}

// position is sample k's position in the window, from −1 to 1
func (f polyFit) position(k int) float64 {
	return float64(k-f.n/2) / float64(f.spacing)
}

// weights maps a window's samples to the fitted polynomial's deriv-th
// derivative, by its position u, at sample at:
// Σ_j q_j(u_k)·q_j^(deriv)(u_at)/‖q_j‖²
func (f polyFit) weights(at, deriv int) (w []float64) {
	w = make([]float64, f.n)
	if deriv > f.polyOrder {
		return // the fit's derivative is 0
	}
	// derivs[m][j] is q_j^(m) at sample at: q_j itself for m = 0,
	// and for m > 0 the recurrence differentiated m times,
	// h[j+1][j]·q_j+1^(m) = u·q_j^(m) + m·q_j^(m−1) − Σ_i≤j h[i][j]·q_i^(m)
	u := f.position(at)
	derivs := make([][]float64, deriv+1)
	derivs[0] = make([]float64, f.polyOrder+1)
	for j := range derivs[0] {
		derivs[0][j] = f.q[j][at]
	}
	for m := 1; m <= deriv; m++ {
		derivs[m] = make([]float64, f.polyOrder+1)
		for j := range f.polyOrder {
			var sum floats.Sum
			sum.Add(float64(u * derivs[m][j]))
			sum.Add(float64(float64(m) * derivs[m-1][j]))
			for i := 0; i <= j; i++ {
				sum.Add(-float64(f.h[i][j] * derivs[m][i]))
			}
			derivs[m][j+1] = sum.Total() / f.h[j+1][j]
		}
	}
	coeffs := make([]float64, f.polyOrder+1)
	for j := deriv; j <= f.polyOrder; j++ {
		coeffs[j] = derivs[deriv][j] / f.norm[j]
	}
	for k := range w {
		var sum floats.Sum
		for j := deriv; j <= f.polyOrder; j++ {
			product := float64(f.q[j][k] * coeffs[j])
			sum.Add(product)
			sum.Add(math.FMA(f.q[j][k], coeffs[j], -product))
		}
		w[k] = sum.Total()
	}
	return
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

// compensatedDot is Σ w_k·x_k,
// each product's rounding error recovered by FMA and summed with compensation,
// so it's as if summed exactly and rounded once,
// except where the terms cancel to far below their size.
// An infinite product has no rounding error to recover, and is summed as it is.
func compensatedDot(w, x []float64) float64 {
	var sum floats.Sum
	for k := range w {
		product := float64(w[k] * x[k])
		sum.Add(product)
		if !math.IsInf(product, 0) {
			sum.Add(math.FMA(w[k], x[k], -product))
		}
	}
	return sum.Total()
}
