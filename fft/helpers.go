package fft

import (
	"fmt"
	"math"
)

// Freq returns the frequency of each of Forward's n bins, for samples d apart:
// 0, 1/(n·d), 2/(n·d), ... up to the highest positive frequency,
// then the negative frequencies, from the most negative up.
// For even n the Nyquist bin is listed as negative.
//
//	f, err := fft.Freq(8, 1.0/8) // [0 1 2 3 -4 -3 -2 -1]
//
// Units:
//   - d is the time between samples, e.g. 1/fs in seconds
//   - f is in cycles per unit of d, so Hz for d in seconds, and cycles per sample for d = 1
//
// Returns:
//   - ErrInvalidLength for a negative n.
//   - ErrInvalidParameter for a d that isn't positive and finite.
//
// Equivalent to:
//   - scipy.fft.fftfreq(n, d)
//   - MATLAB ((0:n-1) - n*((0:n-1) >= ceil(n/2)))/(n*d), as it has no fftfreq
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

// RealFreq returns the frequency of each of ForwardReal's n/2 + 1 bins, for samples d apart:
// 0, 1/(n·d), 2/(n·d), ... up to Nyquist when n is even.
//
//	f, err := fft.RealFreq(len(x), 1/fs) // in Hz
//
// Units:
//   - d is the time between samples, e.g. 1/fs in seconds
//   - f is in cycles per unit of d, so Hz for d in seconds, and cycles per sample for d = 1
//
// Returns:
//   - ErrInvalidLength for a negative n.
//   - ErrInvalidParameter for a d that isn't positive and finite.
//
// Equivalent to:
//   - scipy.fft.rfftfreq(n, d)
//   - MATLAB (0:floor(n/2))/(n*d)
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

// Shift reorders a spectrum, or its frequencies from Freq, so zero frequency is in the centre:
// the negative frequencies first, then the non-negative.
// InverseShift undoes it.
//
//	fft.Shift([]float64{0, 1, 2, 3, -4, -3, -2, -1}) // [-4 -3 -2 -1 0 1 2 3]
//
// Units: same as x
//
// Equivalent to:
//   - scipy.fft.fftshift(x)
//   - MATLAB fftshift(x)
func Shift[T any](x []T) (y []T) {
	y = rotated(x, len(x)-len(x)/2)
	return
}

// InverseShift undoes Shift, putting zero frequency back first.
// For an even length it's the same as Shift, and for an odd length it isn't.
//
// Units: same as x
//
// Equivalent to:
//   - scipy.fft.ifftshift(x)
//   - MATLAB ifftshift(x)
func InverseShift[T any](x []T) (y []T) {
	y = rotated(x, len(x)/2)
	return
}

// FastLength returns the smallest length of at least n whose only prime factors are 2, 3 and 5,
// the lengths tukey's transforms are fastest at.
// Zero-pad to it when any length of at least n will do, e.g. for a correlation:
//
//	y, err := fft.Forward(x, fft.N(fft.FastLength(len(x))))
//
// Units: n and length are numbers of samples
//
// Equivalent to:
//   - scipy.fft.next_fast_len(n, real=True), whose default also allows factors of 7 and 11
//   - MATLAB has none
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
// and ErrInvalidParameter for a sample spacing d that isn't positive and finite.
func checkFreq(n int, d float64) (err error) {
	switch {
	case n < 0:
		err = fmt.Errorf("%w: n %d, must be at least 0", ErrInvalidLength, n)
	case !(d > 0) || math.IsInf(d, 1):
		err = fmt.Errorf("%w: d %g, must be positive and finite", ErrInvalidParameter, d)
	}
	return
}
