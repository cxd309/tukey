package fft

import (
	"fmt"
	"math/cmplx"
	"slices"
)

// ForwardReal returns the discrete Fourier tranform of the real signal x
// keeping on the n/2 + 1 non-negative frequency bins
// a real signal spectrum is Hermitian, X[n-k] = conj(X[k]), so other bins carry no new information
// the DC bin and Nyquist bin when n is even are real: imaginary parts exactly zero
//
// returns ErrInvalidLength if x is empty
// returns ErrInvalidOption for an unknown normalisation
//
// equivalent to scipy.fft.rfft(x)
// MATLAB fft(x) returns all n bins, this is the first n/2 + 1
func ForwardReal(x []float64, opts ...Option) (y []complex128, err error) {
	cfg, err := newConfig(opts, len(x))
	if err != nil {
		return nil, err
	}
	n := cfg.n
	y = slices.Clone(transform(complexified(resized(x, n)))[:n/2+1])

	// exactly real, as they are for any real signal
	// rounding inside the transform can otherwise leave ~1e-16 imaginary parts
	y[0] = complex(real(y[0]), 0)
	if n%2 == 0 {
		y[n/2] = complex(real(y[n/2]), 0)
	}
	scaled(y, cfg.norm.scale(n, true))
	return
}

// InverseReal returns the real signal of length n whose non-negative frequency bins are x
// as ForwardReal produces
// so Inversereal(ForwardReal(x), len(x)) recovers x
// n must be given: n/2 + 1 bins come from both an even and an odd length
// the first n/2 + 1 values of x are used, zero-padded if x is shorter
// the imaginary parts of the DC bin, and the nyquist bin when n is even are ignored
//
// returns ErrInvalidLength if n<1
// returns ErrInvalidOption for an unknown normalisation or if given N
//
// equivalent to scipy.fft.irfft(x, n=n)
// MATLAB ifft(X, 'symmetric') does the same from the full spectrum
func InverseReal(x []complex128, n int, opts ...Option) (y []float64, err error) {
	cfg, err := newConfig(opts, n)
	if err != nil {
		return nil, err
	}
	if cfg.nSet {
		return nil, fmt.Errorf("%w: InverseReal's length is its n argument, not an N option", ErrInvalidOption)
	}
	z := backward(hermitian(x, n))
	scaled(z, cfg.norm.scale(n, false))
	y = make([]float64, n)
	for i := range y {
		y[i] = real(z[i])
	}
	return
}

// hermitian builds the full length-n spectrum of a real signal
// from its non-negative frequency bins:
// X[k] = x[k] and X[n-k] = conj(x[k])
// with x zero-padded or truncated to n/2. + 1 values
// and the DC (for even n) Nyquist bins made real
func hermitian(x []complex128, n int) (full []complex128) {
	half := make([]complex128, n/2+1)
	copy(half, x)

	full = make([]complex128, n)
	full[0] = complex(real(half[0]), 0)
	for k := 1; k < len(half); k++ {
		full[k] = half[k]
		full[n-k] = cmplx.Conj(half[k])
	}
	if n%2 == 0 {
		full[n/2] = complex(real(half[n/2]), 0)
	}
	return
}
