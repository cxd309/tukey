package window

import "math"

// Boxcar returns the n-point rectangular window: all ones
//
// returns ErrInvalidLength if n < 0
//
// equivalent to scipy.signal.windows.boxcar(n),
// equivalent to MATLAB rectwin(n)
func Boxcar(n int, opts ...Option) (w []float64, err error) {
	w, err = newWindow("Boxcar", n, opts, boxcar)
	return
}

// Triang returns the n-point trangular window
// unlike Bartlett, it doesn't reach zero: its ends are 1/n, or 2/(n+1) for odd n
//
// returns ErrInvalidLength if n < 0
//
// equivalent to scipy.signal.windows.trang(n),
// equivalent to MATLAB triang(n)
func Triang(n int, opts ...Option) (w []float64, err error) {
	w, err = newWindow("Triang", n, opts, triang)
	return
}

// Barlett returns the n-point Bartlett window: a triangle, zero at both ends
//
// returns ErrInvalidLength if n < 0
//
// equivalent to scipy.windows.barlett(n),
// equivalent to MATLAB bartlett(n)
func Bartlett(n int, opts ...Option) (w []float64, err error) {
	w, err = newWindow("Barlett", n, opts, bartlett)
	return
}

// Parzen returns the n-point Parzen window:
// a piecewise cubic, the triangle convolved with itself twice
//
// returns ErrInvalidLength if n < 0
//
// equivalent to scipy.signal.windows.parzen(n),
// equivalent to MATLAB parzenwin(n)
func Parzen(n int, opts ...Option) (w []float64, err error) {
	w, err = newWindow("Parzen", n, opts, parzen)
	return
}

// Bohman retuns the n-point Bohman window: a triangle multiplied by a cosine lobe,
// plus a sine term making its first derivative continuous
//
// returns ErrInvalidLength if n < 0
//
// equivalent to scipy.signal.windows.bohman(n),
// equivalent to MATLAB bohmanwin(n)
func Bohman(n int, opts ...Option) (w []float64, err error) {
	w, err = newWindow("Bohman", n, opts, bohman)
	return
}

// Barthann returns the n-point modified Bartlett-Hann window:
// a weighted combination of Bartlett and Hann
//
// returns ErrInvalidLength if n < 0
//
// equivalent to scipy.signal.windows.barthann(n),
// equivalent to MATLAB barthannwindow(n)
func Barthann(n int, opts ...Option) (w []float64, err error) {
	w, err = newWindow("Barthann", n, opts, barthann)
	return
}

// each of these gives a sample k of the first half of a symmetric window of length m,
// which is all shape evaluates

func boxcar(_, _ int) (v float64) {
	v = 1
	return
}

func triang(k, m int) (v float64) {
	if m%2 == 0 {
		v = float64(2*k+1) / float64(m)
	} else {
		v = float64(2*k+2) / float64(m+1)
	}
	return
}

func bartlett(k, m int) (v float64) {
	v = 2 * float64(k) / float64(m-1)
	return
}

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
func bohman(k, m int) (v float64) {
	if k == 0 {
		return 0 // exactly, same as SciPy
	}
	f := 1 - 2*float64(k)/float64(m-1) // distance from the centre, 1 at ends
	v = (1-f)*math.Cos(math.Pi*f) + math.Sin(math.Pi*f)/math.Pi
	return
}

func barthann(k, m int) (v float64) {
	f := math.Abs(float64(k)/float64(m-1) - 0.5)
	v = 0.62 - 0.48*f + 0.38*math.Cos(2*math.Pi*f)
	return
}
