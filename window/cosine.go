package window

import (
	"fmt"
	"math"
)

// GeneralCosine returns the n-point cosine-sum window with coefficients a:
//
//	w(x) = Σ a[j]·cos(j·x), for x evenly spaced from -π to π
//
// It isn't normalised, so its peak is the sum of a, same as SciPy's.
// Hann, Hamming, Blackman, BlackmanHarris, Nuttall and FlatTop are all cosine sums.
//
//	w, err := window.GeneralCosine(64, []float64{1, 1.942604, 1.340318, 0.440811, 0.043097}) // HFT90D
//
// Units:
//   - n is a number of samples
//   - a and w are dimensionless
//
// Returns:
//   - ErrInvalidLength for a negative n.
//   - ErrInvalidParameter for an empty a, or one with a NaN or infinite coefficient.
//   - ErrInvalidOption for an option other than Periodic.
//
// Equivalent to:
//   - scipy.signal.windows.general_cosine(n, a)
//   - MATLAB has none
func GeneralCosine(n int, a []float64, opts ...Option) (w []float64, err error) {
	if len(a) == 0 {
		return nil, fmt.Errorf("%w: GeneralCosine needs at least one coefficient", ErrInvalidParameter)
	}
	w, err = cosineWindow("GeneralCosine", n, a, opts)
	return
}

// GeneralHamming returns the n-point generalised Hamming window:
//
//	w[k] = alpha - (1 - alpha)·cos(2π·k/(n - 1))
//
// alpha 0.5 is Hann, and 0.54 is Hamming.
//
// Units:
//   - n is a number of samples
//   - alpha and w are dimensionless
//
// Returns:
//   - ErrInvalidLength for a negative n.
//   - ErrInvalidParameter for a NaN or infinite alpha.
//   - ErrInvalidOption for an option other than Periodic.
//
// Equivalent to:
//   - scipy.signal.windows.general_hamming(n, alpha)
//   - MATLAB has none
func GeneralHamming(n int, alpha float64, opts ...Option) (w []float64, err error) {
	w, err = cosineWindow("GeneralHamming", n, generalHamming(alpha), opts)
	return
}

// FlatTop returns the n-point flat top window, for measuring a sine's amplitude.
// Its passband is flat to within a fraction of a dB, so a sine between two bins reads its true height.
// It dips below 0, and peaks at 1.000000003, the sum of its published coefficients.
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
//   - scipy.signal.windows.flattop(n)
//   - MATLAB flattopwin(n)
func FlatTop(n int, opts ...Option) (w []float64, err error) {
	w, err = cosineWindow("FlatTop", n, flatTopCoefficients, opts)
	return
}

// Nuttall returns the n-point minimum 4-term Blackman-Harris window, in Nuttall's version.
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
//   - scipy.signal.windows.nuttall(n)
//   - MATLAB nuttallwin(n)
func Nuttall(n int, opts ...Option) (w []float64, err error) {
	w, err = cosineWindow("Nuttall", n, nuttallCoefficients, opts)
	return
}

// BlackmanHarris returns the n-point minimum 4-term Blackman-Harris window,
// whose sidelobes are about 92 dB down.
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
//   - scipy.signal.windows.blackmanharris(n)
//   - MATLAB blackmanharris(n)
func BlackmanHarris(n int, opts ...Option) (w []float64, err error) {
	w, err = cosineWindow("BlackmanHarris", n, blackmanHarrisCoefficients, opts)
	return
}

// Blackman returns the n-point Blackman window:
// three cosine terms, zero at both ends.
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
//   - scipy.signal.windows.blackman(n)
//   - MATLAB blackman(n)
func Blackman(n int, opts ...Option) (w []float64, err error) {
	w, err = cosineWindow("Blackman", n, blackmanCoefficients, opts)
	return
}

// Hamming returns the n-point Hamming window:
// a raised cosine on a pedestal of 0.08.
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
//   - scipy.signal.windows.hamming(n)
//   - MATLAB hamming(n)
func Hamming(n int, opts ...Option) (w []float64, err error) {
	w, err = cosineWindow("Hamming", n, generalHamming(0.54), opts)
	return
}

// Hann returns the n-point Hann window:
// a raised cosine, zero at both ends, or only at the start when periodic.
//
//	w, err := window.Hann(64, window.Periodic()) // for spectral analysis
//
// MATLAB's hanning(n) is Hann without its zero ends, Hann(n+2) less its first and last sample:
//
//	w, err := window.Hann(n + 2)
//	w = w[1 : n+1]
//
// MATLAB's hanning(n, 'periodic') is Hann(n, Periodic()).
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
//   - scipy.signal.windows.hann(n), and numpy.hanning(n)
//   - MATLAB hann(n)
func Hann(n int, opts ...Option) (w []float64, err error) {
	w, err = cosineWindow("Hann", n, generalHamming(0.5), opts)
	return
}

// The coefficients of the fixed cosine-sum windows, exactly as SciPy's.
// FlatTop's sum to 1.000000003, as published, and MATLAB's flattopwin uses the same.
var (
	blackmanCoefficients       = []float64{0.42, 0.50, 0.08}
	blackmanHarrisCoefficients = []float64{0.35875, 0.48829, 0.14128, 0.01168}
	nuttallCoefficients        = []float64{0.3635819, 0.4891775, 0.1365995, 0.0106411}
	flatTopCoefficients        = []float64{0.21557895, 0.41663158, 0.277263158, 0.083578947, 0.006947368}
)

// cosineSum gives sample k of the symmetric cosine-sum window of length m:
// the sum of a[j]·cos(j·x), for x evenly spaced from -π to π.
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

// cosineWindow builds the n-point cosine-sum window with coefficients a,
// returning ErrInvalidParameter for a NaN or infinite coefficient.
func cosineWindow(name string, n int, a []float64, opts []Option) (w []float64, err error) {
	for _, v := range a {
		if math.IsNaN(v) || math.IsInf(v, 0) {
			return nil, fmt.Errorf("%w: %s's coefficients must be finite, got %v", ErrInvalidParameter, name, a)
		}
	}
	w, err = newWindow(name, n, opts, cosineSum(a))
	return
}

// generalHamming is the cosine-sum coefficients of the generalised Hamming window.
func generalHamming(alpha float64) (a []float64) {
	a = []float64{alpha, 1 - alpha}
	return
}
