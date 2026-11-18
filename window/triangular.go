package window

import "math"

// Boxcar returns the n-point rectangular window, all ones.
//
// Units:
//   - n is a number of samples
//   - w is dimensionless
//
// Returns:
//   - ErrInvalidLength for a negative n.
//   - ErrInvalidOption for an option other than Periodic.
//
// Equivalent to:
//   - scipy.signal.windows.boxcar(n)
//   - MATLAB rectwin(n), which has no periodic form
func Boxcar(n int, opts ...Option) (w []float64, err error) {
	w, err = newWindow("Boxcar", n, opts, boxcar)
	return
}

// Triang returns the n-point triangular window.
// Unlike Bartlett it doesn't reach 0, as its ends are 1/n, or 2/(n+1) for odd n.
//
// Units:
//   - n is a number of samples
//   - w is dimensionless
//
// Returns:
//   - ErrInvalidLength for a negative n.
//   - ErrInvalidOption for an option other than Periodic.
//
// Equivalent to:
//   - scipy.signal.windows.triang(n)
//   - MATLAB triang(n), which has no periodic form
func Triang(n int, opts ...Option) (w []float64, err error) {
	w, err = newWindow("Triang", n, opts, triang)
	return
}

// Bartlett returns the n-point Bartlett window, a triangle that's 0 at both ends.
//
// Units:
//   - n is a number of samples
//   - w is dimensionless
//
// Returns:
//   - ErrInvalidLength for a negative n.
//   - ErrInvalidOption for an option other than Periodic.
//
// Equivalent to:
//   - scipy.signal.windows.bartlett(n)
//   - MATLAB bartlett(n), which has no periodic form
func Bartlett(n int, opts ...Option) (w []float64, err error) {
	w, err = newWindow("Bartlett", n, opts, bartlett)
	return
}

// Parzen returns the n-point Parzen window:
// a piecewise cubic, a triangle convolved with itself.
//
// Units:
//   - n is a number of samples
//   - w is dimensionless
//
// Returns:
//   - ErrInvalidLength for a negative n.
//   - ErrInvalidOption for an option other than Periodic.
//
// Equivalent to:
//   - scipy.signal.windows.parzen(n)
//   - MATLAB parzenwin(n), which has no periodic form
func Parzen(n int, opts ...Option) (w []float64, err error) {
	w, err = newWindow("Parzen", n, opts, parzen)
	return
}

// Bohman returns the n-point Bohman window:
// a triangle multiplied by a cosine lobe, plus a sine term that makes its first derivative continuous.
//
// Units:
//   - n is a number of samples
//   - w is dimensionless
//
// Returns:
//   - ErrInvalidLength for a negative n.
//   - ErrInvalidOption for an option other than Periodic.
//
// Equivalent to:
//   - scipy.signal.windows.bohman(n)
//   - MATLAB bohmanwin(n), which has no periodic form
func Bohman(n int, opts ...Option) (w []float64, err error) {
	w, err = newWindow("Bohman", n, opts, bohman)
	return
}

// Barthann returns the n-point modified Bartlett-Hann window:
// a weighted combination of Bartlett and Hann.
//
// Units:
//   - n is a number of samples
//   - w is dimensionless
//
// Returns:
//   - ErrInvalidLength for a negative n.
//   - ErrInvalidOption for an option other than Periodic.
//
// Equivalent to:
//   - scipy.signal.windows.barthann(n)
//   - MATLAB barthannwin(n), which has no periodic form
func Barthann(n int, opts ...Option) (w []float64, err error) {
	w, err = newWindow("Barthann", n, opts, barthann)
	return
}

// boxcar is shape's value for the rectangular window, 1 everywhere.
func boxcar(_, _ int) (v float64) {
	v = 1
	return
}

// triang is shape's value for the triangular window, which is nonzero at its ends.
func triang(k, m int) (v float64) {
	if m%2 == 0 {
		v = float64(2*k+1) / float64(m)
	} else {
		v = float64(2*k+2) / float64(m+1)
	}
	return
}

// bartlett is shape's value for the Bartlett window, which is 0 at its ends.
func bartlett(k, m int) (v float64) {
	v = 2 * float64(k) / float64(m-1)
	return
}

// parzen is shape's value for the Parzen window, a piecewise cubic.
func parzen(k, m int) (v float64) {
	fromCentre := float64(m-1)/2 - float64(k)
	a := fromCentre / (float64(m) / 2)
	if fromCentre <= float64(m-1)/4 {
		v = 1 - 6*math.Pow(a, 2) + 6*math.Pow(a, 3)
	} else {
		v = 2 * math.Pow(1-a, 3)
	}
	return
}

// bohman is shape's value for the Bohman window.
func bohman(k, m int) (v float64) {
	if k == 0 {
		return 0 // exactly, same as SciPy
	}
	f := 1 - 2*float64(k)/float64(m-1) // distance from the centre, 1 at the ends
	v = (1-f)*math.Cos(math.Pi*f) + math.Sin(math.Pi*f)/math.Pi
	return
}

// barthann is shape's value for the Bartlett-Hann window.
func barthann(k, m int) (v float64) {
	f := math.Abs(float64(k)/float64(m-1) - 0.5)
	v = 0.62 - 0.48*f + 0.38*math.Cos(2*math.Pi*f)
	return
}
