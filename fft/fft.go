package fft

import (
	"math/cmplx"
)

// Forward returns the discrete Fourier transform of x, n = len(x) bins:
//
//	y[k] = Σ x[j]·e^(-2πi·jk/n), for k = 0, 1, ..., n-1
//
// scaled as the Norm options say, unscaled by default.
// Bin k is at frequency k/n cycles per sample, and Freq gives each bin's frequency.
//
//	y, err := fft.Forward(x)
//
// Units: same as x
//
// Returns:
//   - ErrSignalTooShort for an empty x without N.
//   - ErrInvalidLength for an N below 1.
//
// Equivalent to:
//   - scipy.fft.fft(x)
//   - MATLAB fft(x)
func Forward(x []complex128, opts ...Option) (y []complex128, err error) {
	cfg, err := newConfig(opts, len(x))
	if err != nil {
		return nil, err
	}
	y = transform(resized(x, cfg.n))
	scaled(y, cfg.norm.scale(cfg.n, true))
	return
}

// Inverse returns the inverse discrete Fourier transform of x, n = len(x) samples:
//
//	y[j] = (1/n)·Σ x[k]·e^(+2πi·jk/n), for j = 0, 1, ..., n-1
//
// so Inverse(Forward(x)) recovers x.
// The Norm options move the 1/n.
//
// Units: same as x
//
// Returns:
//   - ErrSignalTooShort for an empty x without N.
//   - ErrInvalidLength for an N below 1.
//
// Equivalent to:
//   - scipy.fft.ifft(x)
//   - MATLAB ifft(x)
func Inverse(x []complex128, opts ...Option) (y []complex128, err error) {
	cfg, err := newConfig(opts, len(x))
	if err != nil {
		return nil, err
	}
	y = backward(resized(x, cfg.n))
	scaled(y, cfg.norm.scale(cfg.n, false))
	return
}

// backward is the unscaled inverse transform:
// the forward transform with conjugated twiddles, which is conj(transform(conj(x))).
func backward(x []complex128) (y []complex128) {
	y = conjugated(transform(conjugated(x)))
	return
}

// conjugated returns a new slice holding the complex conjugate of each value in x.
func conjugated(x []complex128) (c []complex128) {
	c = make([]complex128, len(x))
	for i, v := range x {
		c[i] = cmplx.Conj(v)
	}
	return
}

// transform is the unscaled forward transform, by the fastest algorithm for x's length.
func transform(x []complex128) (y []complex128) {
	n := len(x)
	if isPowerOfTwo(n) {
		y = radix2(x)
	} else {
		y = mixedRadix(x)
	}
	return
}

// complexified returns x as complex values with zero imaginary parts.
func complexified(x []float64) (c []complex128) {
	c = make([]complex128, len(x))
	for i, v := range x {
		c[i] = complex(v, 0)
	}
	return
}
