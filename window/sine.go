package window

import "math"

// Cosine returns the n-point cosine (sine) window:
// half a peroid of a sine, not reaching zero at its ends
//
// returns ErrInvalidLength if n < 0
//
// equivalent to scipy.signal.windows.cosine(n),
// no MATLAB equivalent
func Cosine(n int, opts ...Option) (w []float64, err error) {
	w, err = newWindow("Cosine", n, opts, cosine)
	return
}

// Lanczos returns the n-point Lanczos (sinc) window:
// the central lobe of sin(pi x)/(pi x), for x from -1 to 1
//
// returns ErrInvalidLength if n < 0
//
// equivalent to scipy.signal.windows.lanczos(n),
// no MATLAB equivalent
func Lanczos(n int, opts ...Option) (w []float64, err error) {
	w, err = newWindow("Lanczos", n, opts, lanczos)
	return
}

func cosine(k, m int) (v float64) {
	v = math.Sin(math.Pi / float64(m) * (float64(k) + 0.5))
	return
}

func lanczos(k, m int) (v float64) {
	x := 2*float64(k)/float64(m-1) - 1
	if x == 0 {
		return 1 // sinc(0), the centre of an odd window
	}
	v = math.Sin(math.Pi*x) / (math.Pi * x)
	return
}
