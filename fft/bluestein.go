package fft

import "math/cmplx"

// bluestein computes the DFT of any length n as a convolution,
// which fast transforms of a length of at least 2n - 1 compute.
// Writing jk = (j² + k² - (k - j)²)/2 turns the DFT into:
//
//	X[k] = c[k]·Σ (x[j]·c[j])·conj(c[k - j]), where c[j] = e^(-πi·j²/n)
//
// The sum is the convolution of x·c with conj(c), done by multiplying spectra.
func bluestein(x []complex128) (y []complex128) {
	n := len(x)
	chirp := make([]complex128, n)
	for j := range chirp {
		chirp[j] = twiddle(j*j%(2*n), 2*n) // e^(-πi·j²/n), reduced exactly
	}

	size := paddedSize(n)

	// x·c, zero-padded
	a := make([]complex128, size)
	for j := range n {
		a[j] = x[j] * chirp[j]
	}
	// conj(c) laid out circularly: c is even in j, so index -j wraps to size - j
	b := make([]complex128, size)
	b[0] = cmplx.Conj(chirp[0])
	for j := 1; j < n; j++ {
		b[j] = cmplx.Conj(chirp[j])
		b[size-j] = cmplx.Conj(chirp[j])
	}

	// the convolution theorem: multiply spectra, then transform back
	spectrum, kernel := transform(a), transform(b)
	for i := range spectrum {
		spectrum[i] *= kernel[i]
	}
	convolution := conjugated(transform(conjugated(spectrum))) // inverse, unscaled

	scale := complex(1/float64(size), 0)
	y = make([]complex128, n)
	for k := range y {
		y[k] = chirp[k] * convolution[k] * scale
	}
	return
}

// paddedSize is the convolution length Bluestein uses for an n-point transform:
// at least 2n - 1, and either the next power of two,
// or the next length with no prime factors above 5 when that's short enough to make up for its cost.
// Mixed radix costs about 1.5 times radix-2 per element, 2.0 to 2.5 against 1.5 ns per n·log2(n).
func paddedSize(n int) (size int) {
	target := 2*n - 1
	powerOfTwo := 1
	for powerOfTwo < target {
		powerOfTwo *= 2
	}
	smooth := smoothSize(target)
	if 1.5*float64(smooth) < float64(powerOfTwo) {
		return smooth
	}
	return powerOfTwo
}

// smoothSize returns the smallest length of at least target whose only prime factors are 2, 3 and 5,
// the factors with dedicated butterflies.
// It tries every product of powers of 3 and 5 below the best so far, doubling each up to target.
func smoothSize(target int) (size int) {
	size = 1
	for size < target {
		size *= 2 // the power of two is always a candidate
	}
	for f5 := 1; f5 < size; f5 *= 5 {
		for f35 := f5; f35 < size; f35 *= 3 {
			candidate := f35
			for candidate < target {
				candidate *= 2
			}
			size = min(size, candidate)
		}
	}
	return
}
