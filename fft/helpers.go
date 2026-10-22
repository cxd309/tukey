package fft

import "slices"

// Freq returns the frequency of each of Forward's n bins for samples d apart
// 0, 1/(n d), 2/(n d), ... up to the highest positive frequency
// then the negative frequencies from the most negative frequency up
// With d = 1/fs the frequencies are all in Hz
// For even n, Nyquist bin is listed as negative
//
// equivalent to scipy.fft.fftfreq(n, d)
func Freq(n int, d float64) (f []float64) {
	f = make([]float64, max(n, 0))
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
// equivalent to scipy.fft.rfftfreq(n, d)
func RealFreq(n int, d float64) (f []float64) {
	f = make([]float64, max(n/2+1, 0))
	if n < 1 {
		f = f[:0]
	}
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
	cut := len(x) - (len(x) / 2)
	y = slices.Concat(x[cut:], x[:cut])
	return
}

// InverseShift undoes Shift, putting zero frequencies back first, for odd lengths it differs from Shift
//
// equivalent to scipy.fft.ifftshift and MATLAB ifftshift
func InverseShift[T any](x []T) (y []T) {
	cut := len(x) / 2
	y = slices.Concat(x[cut:], x[:cut])
	return
}

// FastLength returns the smallest length of at least n whose only prime
// factors are 2, 3, 5, the lengths tukey's FFTs transform fastest.
// Zero-pad to it when any length of at least n will do, as for a correlation.
func FastLength(n int) (length int) {
	length = smoothSize(n)
	return
}
