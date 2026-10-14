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
	signal := resized(x, n)
	if n%2 == 0 {
		y = forwardRealPacked(signal)
	} else {
		y = slices.Clone(transform(complexified(signal))[:n/2+1])
		y[0] = complex(real(y[0]), 0) // exactly real, as for any real signal
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
	if n%2 == 0 {
		y = inverseRealPacked(x, n)
	} else {
		z := backward(hermitian(x, n))
		y = make([]float64, n)
		for i := range y {
			y[i] = real(z[i])
		}
	}
	if s := cfg.norm.scale(n, false); s != 1 {
		for i := range y {
			y[i] *= s
		}
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

// forwardRealPacked is ForwardReal's transform for even n, at half the cost of the complex one
// with h = n/2, it packs pairs of samples into one complex value,
// z[j] = x[2j] + i·x[2j+1]
// so one h-point FFT gives Z = E + i·O, where E and O are the transforms of the even and odd samples
// Each is a real sequence's spectrum, so Hermitian, which separates them:
//
//	E[k] = (Z[k] + conj(Z[h-k]))/2,   O[k] = -i·(Z[k] - conj(Z[h-k]))/2
//
// then a radix-2 butterfly with w = e^(-2πi·k/n) combines them, giving bins k
// and h-k from the same pair:
//
//	X[k] = E[k] + w·O[k],   X[h-k] = conj(E[k] - w·O[k])
func forwardRealPacked(x []float64) (y []complex128) {
	h := len(x) / 2
	z := make([]complex128, h)
	for j := range z {
		z[j] = complex(x[2*j], x[2*j+1])
	}
	spectrum := transform(z)

	y = make([]complex128, h+1)
	y[0] = complex(real(spectrum[0])+imag(spectrum[0]), 0) // DC: E[0] + O[0], both real
	y[h] = complex(real(spectrum[0])-imag(spectrum[0]), 0) // Nyquist: E[0] - O[0]
	for k := 1; k <= h/2; k++ {
		a, b := spectrum[k], cmplx.Conj(spectrum[h-k])
		even := scaledBy(0.5, a+b)
		odd := minusI(scaledBy(0.5, a-b))
		rotated := twiddle(k, len(x)) * odd
		y[k] = even + rotated
		y[h-k] = cmplx.Conj(even - rotated)
	}
	return
}

// inverseRealPacked is InverseReal's unscaled transform for even n, returning n·x:
// forwardRealPacked's steps in reverse
// It rebuilds Z = 2E + 2i·O from the bins, without the halves, so the unscaled h-point inverse gives exactly n·x,
// matching the 1/n the normalisation then applies:
//
//	Z[k] = A + i·conj(w)·B,   Z[h-k] = conj(A - i·conj(w)·B)
//	where A = X[k] + conj(X[h-k]), B = X[k] - conj(X[h-k])
//
// then unpacks the even samples from the real parts and the odd from the imaginary.
// Only the real parts of the DC and Nyquist bins are used.
func inverseRealPacked(x []complex128, n int) (y []float64) {
	h := n / 2
	bins := make([]complex128, h+1) // zero-padded or truncated to h+1, as SciPy
	copy(bins, x)
	dc, nyquist := real(bins[0]), real(bins[h])

	spectrum := make([]complex128, h)
	spectrum[0] = complex(dc+nyquist, dc-nyquist)
	for k := 1; k <= h/2; k++ {
		a := bins[k] + cmplx.Conj(bins[h-k])
		b := bins[k] - cmplx.Conj(bins[h-k])
		rotated := complex(0, 1) * cmplx.Conj(twiddle(k, n)) * b
		spectrum[k] = a + rotated
		spectrum[h-k] = cmplx.Conj(a - rotated)
	}

	z := backward(spectrum)
	y = make([]float64, n)
	for j, v := range z {
		y[2*j], y[2*j+1] = real(v), imag(v)
	}
	return
}
