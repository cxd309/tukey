package fft

import "math/cmplx"

// bluestein computes the DFT of any length n as a convolution
// which power of two FFTs of length L>= 2n-1 compute fast
// From jk = (j^2 + k^2 - (k-j)^2) / 2
//
//	X[k] = c[k] sum_j (x[j] c[j]) conj(c[k-j]), c[j] = e^(-pi i j^2 /n)
//
// the sum is the convolution of xc with conj(c) done by multiplying spectra
func bluestein(x []complex128) (y []complex128) {
	n := len(x)
	chirp := make([]complex128, n)
	for j := range chirp {
		chirp[j] = twiddle(j*j%(2*n), 2*n) // e^(-pi i j^2 /n) reduced exactly
	}

	size := 1
	for size < 2*n-1 {
		size *= 2
	}

	// xc zero-padded
	a := make([]complex128, size)
	for j := range n {
		a[j] = x[j] * chirp[j]
	}
	// conj(c) laid out circularly c is even in j so index -j wraps to size-j
	b := make([]complex128, size)
	b[0] = cmplx.Conj(chirp[0])
	for j := 1; j < n; j++ {
		b[j] = cmplx.Conj(chirp[j])
		b[size-j] = cmplx.Conj(chirp[j])
	}

	// the convolution theorem: multiply spectra, then transform back
	spectrum, kernel := radix2(a), radix2(b)
	for i := range spectrum {
		spectrum[i] *= kernel[i]
	}
	convolution := conjugated(radix2(conjugated(spectrum))) // inverse, unscaled

	scale := complex(1/float64(size), 0)
	y = make([]complex128, n)
	for k := range y {
		y[k] = chirp[k] * convolution[k] * scale
	}
	return
}
