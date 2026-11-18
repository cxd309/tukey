package window

import "math"

// Cosine returns the n-point cosine, or sine, window:
// half a period of a sine, not reaching 0 at its ends.
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
//   - scipy.signal.windows.cosine(n)
//   - MATLAB has none
func Cosine(n int, opts ...Option) (w []float64, err error) {
	w, err = newWindow("Cosine", n, opts, cosine)
	return
}

// Lanczos returns the n-point Lanczos, or sinc, window:
// the central lobe of sin(π·x)/(π·x), for x from -1 to 1.
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
//   - scipy.signal.windows.lanczos(n)
//   - MATLAB has none
func Lanczos(n int, opts ...Option) (w []float64, err error) {
	w, err = newWindow("Lanczos", n, opts, lanczos)
	return
}

// cosine is shape's value for the cosine window, half a period of a sine.
func cosine(k, m int) (v float64) {
	v = math.Sin(math.Pi / float64(m) * (float64(k) + 0.5))
	return
}

// lanczos is shape's value for the Lanczos window, the central lobe of a sinc.
func lanczos(k, m int) (v float64) {
	x := 2*float64(k)/float64(m-1) - 1
	if x == 0 {
		return 1 // sinc(0), the centre of an odd window
	}
	v = math.Sin(math.Pi*x) / (math.Pi * x)
	return
}
