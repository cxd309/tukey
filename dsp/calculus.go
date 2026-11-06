package dsp

import (
	"fmt"
	"math"

	"github.com/cxd309/tukey/internal/floats"
)

// Diff returns the difference between each sample and the next,
// d[i] = x[i+1] - x[i], with a length one short than x.
// Fewer than 2 samples have no differences, so return an empty slice.
//
// To call a higher order use Diff(Diff(x)) etc.
//
// Equivalent to numpy.diff(x).
// Equivalent to MATLAB diff(x).
func Diff(x []float64) (d []float64) {
	d = make([]float64, max(len(x)-1, 0))
	for i := range d {
		d[i] = x[i+1] - x[i]
	}
	return
}

// Grad returns the derivative of x estimated at each sample, sampled every dx:
// central differences inside, (x[i+1] - x[i-1]) / (2·dx), and one-sided differences at the ends.
// Where Diff gives the slope between samples, Grad gives it at each one, the same length as x.
//
// Returns
//   - ErrInvalidLength for fewer than 2 samples
//   - ErrInvalidParameter for a dx that isn't positive and finite
//
// Equivalent to:
//   - numpy.gradient(x, dx)
//   - MATLAB gradient(x, dx)
func Grad(x []float64, dx float64) (g []float64, err error) {
	if len(x) < 2 {
		return nil, fmt.Errorf("%w: Gerad of %d samples, needs at least 2", ErrInvalidLength, len(x))
	}
	if err = checkSpacing(dx); err != nil {
		return nil, err
	}
	n := len(x)
	g = make([]float64, n)
	g[0] = (x[1] - x[0]) / dx
	for i := 1; i < n-1; i++ {
		g[i] = (x[i+1] - x[i-1]) / (2 * dx)
	}
	g[n-1] = (x[n-1] - x[n-2]) / dx
	return
}

// CumSum returns the running total of x, s[i] = x[0] + ... + x[i], the same length as x.
// It's summed with compensation, so its error doesn't grow with the signal's length as a plain running sum's does.
//
// Equivalent to:
//   - numpy.cumsum(x)
//   - MATLAB cumsum(x)
//
// both of which sum plainly
func CumSum(x []float64) (s []float64) {
	s = make([]float64, len(x))
	var sum floats.Sum
	for i, v := range x {
		sum.Add(v)
		s[i] = sum.Total()
	}
	return
}

// Trapz returns the trapezoidal integral of y sampled every dx,
// over the whole signal: dx·(y[0]/2 + y[1] + ... + y[n-1]/2), summed with compensation.
// Fewer than 2 samples have no interval to integrate over, so give 0.
//
// Returns ErrInvalidParameter for a dx that isn't positive and finite.
//
// Equivalent to:
//   - MATLAB trapz(dx, y)
//   - scipy.integrate.trapezoid(y, dx=0)
//   - numpy.trapezoid (renamed from trapz in NumPy 2)
func Trapz(y []float64, dx float64) (integral float64, err error) {
	if err = checkSpacing(dx); err != nil {
		return
	}
	if len(y) < 2 {
		return
	}
	var sum floats.Sum
	sum.Add(y[0] / 2)
	for _, v := range y[1 : len(y)-1] {
		sum.Add(v)
	}
	sum.Add(y[len(y)-1] / 2)
	integral = dx * sum.Total()
	return
}

// CumTrapz returns the running trapezoidal integral of y sampled every dx:
// the integral up to each sample, starting at 0, the same length as y, so it lines up with y and its time axis.
// Summed with compensation, same as Trapz.
//
// Returns ErrInvalidParameter for a dx that isn't positive and finite
//
// Equivalent to:
//   - MATLAB cumtrapz(dx, y)
//   - scipy.integrate.cumulative_trapezoid(y, dx=dx, initial=0)
func CumTrapz(y []float64, dx float64) (integral []float64, err error) {
	if err = checkSpacing(dx); err != nil {
		return nil, err
	}
	integral = make([]float64, len(y))
	var sum floats.Sum // y[0]/2 + y[1] + ... + y[i-1] + y[i]/2
	for i := 1; i < len(y); i++ {
		sum.Add(y[i-1] / 2)
		sum.Add(y[i] / 2)
		integral[i] = dx * sum.Total()
	}
	return
}

// checkSpacing returns ErrInvalidParameter for a sample spacing dx that isn't positive and finite
func checkSpacing(dx float64) (err error) {
	if !(dx > 0) {
		err = fmt.Errorf("%w: dx %g, must be positive", ErrInvalidParameter, dx)
	}
	if math.IsInf(dx, 1) {
		err = fmt.Errorf("%w: dx %g, must be finite", ErrInvalidParameter, dx)
	}
	return
}
