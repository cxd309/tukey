package fft

import (
	"fmt"
	"math"
)

// Freq returns the frequency of each of Forward's n bins for samples d apart
// 0, 1/(n d), 2/(n d), ... up to the highest positive frequency
// then the negative frequencies from the most negative frequency up
// With d = 1/fs the frequencies are all in Hz
// For even n, Nyquist bin is listed as negative
//
// Returns:
//   - ErrInvalidLength for a negative n.
//   - ErrInvalidParameter for a d that isn't positive and finite.
//
// equivalent to scipy.fft.fftfreq(n, d)
func Freq(n int, d float64) (f []float64, err error) {
	if err = checkFreq(n, d); err != nil {
		return nil, err
	}
	f = make([]float64, n)
	step := 1 / (float64(n) * d)
	positive := (n-1)/2 + 1 // bins 0...positive-1 are non-negative
	for i := range f {
		k := i
		if i >= positive {
			k = i - n
		}
		f[i] = float64(k) * step
	}
	return
}

// Realfreq returns the frequency of each of ForwardReal's n/2 + 1 bins for samples d apart
// 0, 1/(n d), 2/(n d)... up to Nyquist when n is even
//
// Returns:
//   - ErrInvalidLength for a negative n.
//   - ErrInvalidParameter for a d that isn't positive and finite.
//
// equivalent to scipy.fft.rfftfreq(n, d)
func RealFreq(n int, d float64) (f []float64, err error) {
	if err = checkFreq(n, d); err != nil {
		return nil, err
	}
	if n == 0 {
		return []float64{}, nil
	}
	f = make([]float64, n/2+1)
	step := 1 / (float64(n) * d)
	for i := range f {
		f[i] = float64(i) * step
	}
	return
}

// Shift reorders a spectrum (or its frequencies, from Freq) so zero frequency is in the centre
// the negative freqencies first, the non-negative
// undone by InverseShift
//
// equivalent to scipy.fft.fftshift and MATLAB fftshift
func Shift[T any](x []T) (y []T) {
	y = rotated(x, len(x)-len(x)/2)
	return
}

// InverseShift undoes Shift, putting zero frequencies back first, for odd lengths it differs from Shift
//
// equivalent to scipy.fft.ifftshift and MATLAB ifftshift
func InverseShift[T any](x []T) (y []T) {
	y = rotated(x, len(x)/2)
	return
}

// FastLength returns the smallest length of at least n whose only prime
// factors are 2, 3, 5, the lengths tukey's FFTs transform fastest.
// Zero-pad to it when any length of at least n will do, as for a correlation.
func FastLength(n int) (length int) {
	length = smoothSize(n)
	return
}

// rotated is a new slice holding x from cut to its end, then from its start to cut.
// It's never nil, so an empty x gives an empty slice.
func rotated[T any](x []T, cut int) (y []T) {
	y = make([]T, len(x))
	n := copy(y, x[cut:])
	copy(y[n:], x[:cut])
	return
}

// checkFreq returns ErrInvalidLength for a negative n,
// and ErrInvalidParameter for a sample spacing d that isn't positive and finite
func checkFreq(n int, d float64) (err error) {
	switch {
	case n < 0:
		err = fmt.Errorf("%w: n %d, must be at least 0", ErrInvalidLength, n)
	case !(d > 0) || math.IsInf(d, 1):
		err = fmt.Errorf("%w: d %g, must be positive and finite", ErrInvalidParameter, d)
	}
	return
}
