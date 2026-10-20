package window

import (
	"fmt"
	"math"
)

// coefficients of the fixed cosine-sum windows, exactly as SciPy's
// (flat top's sum to 1.000000003, as published. MATLAB's flattopwin uses the same)
var (
	blackmanCoefficients       = []float64{0.42, 0.50, 0.08}
	blackmanHarrisCoefficients = []float64{0.35875, 0.48829, 0.14128, 0.01168}
	nuttallCoefficients        = []float64{0.3635819, 0.4891775, 0.1365995, 0.0106411}
	flatTopCoefficients        = []float64{0.21557895, 0.41663158, 0.277263158, 0.083578947, 0.006947368}
)

// cosineSum gives sample k of the symmetric cosine-sum window of length m
// the the sum of a[j] cos(j x), for x evenly spaced from -pi to pi
func cosineSum(a []float64) (value func(k, m int) float64) {
	value = func(k, m int) (v float64) {
		x := 2*math.Pi*float64(k)/float64(m-1) - math.Pi
		for j, aj := range a {
			v += aj * math.Cos(float64(j)*x)
		}
		return
	}
	return
}

// cosineWindow builds the n-point cosine-sum window with coefficients a
func cosineWindow(name string, n int, a []float64, opts []Option) (w []float64, err error) {
	w, err = newWindow(name, n, opts, cosineSum(a))
	return
}

// generalHamming is the cosine-sum coefficients of the generalised Hamming window
func generalHamming(alpha float64) (a []float64) {
	a = []float64{alpha, 1 - alpha}
	return
}

// GeneralCosine returns the n-point cosine-sum window with coefficients a:
// the sum of a[j] cos(j x), for x evenly spaced from -pi to pi
// it isn't normalised: its peak is the sum of a, as SciPy's
//
// returns ErrInvalidLength if n < 0, and ErrInvalidParameter if a is empty
//
// equivalent to scipy.signal.windows.general_cosine(n, a)
// no MATLAB equivalent
func GeneralCosine(n int, a []float64, opts ...Option) (w []float64, err error) {
	if len(a) == 0 {
		return nil, fmt.Errorf("%w: GeneralCosine needs at least one coefficient", ErrInvalidParameter)
	}
	w, err = cosineWindow("GeneralCosine", n, a, opts)
	return
}

// Hanning returns MATLAB's n-point Hanning window:
// symmetric, it's the Hann window without the zero endpoints,  Hann(n+2) less its first and last sample.
// periodic, it's Hann's periodic window, keeping it's first zero
//
// returns ErrInvalidLength if n<0
//
// equivalent to MATLAB's hanning(n), no SciPy equivalent
func Hanning(n int, opts ...Option) (w []float64, err error) {
	cfg, err := newConfig("Hanning", opts, optionSet{optionPeriodic: true})
	if err != nil {
		return nil, err
	}
	hann := cosineSum(generalHamming(0.5))
	if cfg.periodic {
		w, err = shape(n, true, hann)
		return
	}
	w, err = shape(n, false, func(k, m int) float64 { return hann(k+1, m+2) })
	return
}

// GeneralHamming returns the n-point generalised Hamming window
// alpha - (1 - alpha) cos(2 pi k / (n-1))
// alpha 0.5 is Hann, 0.54 is Hamming
//
// returns ErrInvalidLength if n < 0
//
// equivalnet to scipy.signal.windows.general_hamming(n, alpha)
// no MATLAB equivalent
func GeneralHamming(n int, alpha float64, opts ...Option) (w []float64, err error) {
	w, err = cosineWindow("GeneralHamming", n, generalHamming(alpha), opts)
	return
}

// FlatTop returns the n-point float top window, for measuring amp. of sins
// its passband is flat to within a fraction of a dB
// it dips below zero and peaks at 1.000000003
//
// returns ErrInvalidLength if n < 0
//
// equivalent to scipy.signal.windows.flattop(n) and MATLAB's flattopwin(n)
func FlatTop(n int, opts ...Option) (w []float64, err error) {
	w, err = cosineWindow("FlatTop", n, flatTopCoefficients, opts)
	return
}

// Nuttall returns the n-point minimum 4-term Blackman–Harris window, Nuttall's version
//
// returns ErrInvalidLength if n < 0
//
// equivalent to scipy.signal.windows.nuttall(n)
// equivalent to MATLAB's nuttallwin(n)
func Nuttall(n int, opts ...Option) (w []float64, err error) {
	w, err = cosineWindow("Nuttall", n, nuttallCoefficients, opts)
	return
}

// BlackmanHarris returns the n-point minimum 4-term Blackman–Harris window
//
// returns ErrInvalidLength if n < 0
//
// equivalent to scipy.signal.windows.blackmanharris(n)
// equivalent to MATLAB's blackmanharris(n)
func BlackmanHarris(n int, opts ...Option) (w []float64, err error) {
	w, err = cosineWindow("BlackmanHarris", n, blackmanHarrisCoefficients, opts)
	return
}

// Blackman returns the n-point Blackman window: three cosine terms, zero at both ends
//
// returns ErrInvalidLength if n < 0
//
// equivalent to scipy.signal.windows.blackman(n)
// equivalent to MATLAB's blackman(n)
func Blackman(n int, opts ...Option) (w []float64, err error) {
	w, err = cosineWindow("Blackman", n, blackmanCoefficients, opts)
	return
}

// Hamming returns the n-point Hamming window: a raised cosine on a pedestal of 0.08
//
// returns ErrInvalidLength if n < 0
//
// equivalent to scipy.signal.windows.hamming(n)
// equivalent to MATLAB's hamming(n)
func Hamming(n int, opts ...Option) (w []float64, err error) {
	w, err = cosineWindow("Hamming", n, generalHamming(0.54), opts)
	return
}

// Hann returns the n-point Hann window:
// a raised cosine, zero at both ends, or only at the start when periodic
//
// returns ErrInvalidLength if n < 0
//
// equivalent to scipy.signal.windows.hann(n)
// equivalnt to MATLAB's hann(n)
func Hann(n int, opts ...Option) (w []float64, err error) {
	w, err = cosineWindow("Hann", n, generalHamming(0.5), opts)
	return
}
