package fft

import (
	"fmt"
	"math/cmplx"
)

// Forward returns the discrete Fourier transform of x
// y[k] = sum x[j] e^(-2 pi i j k/n) for k = 0 ... n-1
//
// returns ErrInvalidLength if x is empty
//
// equivalent to scipy.fft.fft(x) and MATLAB's fft(x)
func Forward(x []complex128) (y []complex128, err error) {
	if err = checkLength(len(x)); err != nil {
		return nil, err
	}
	y = transform(x)
	return
}

// Inverse returns the inverse discrete Fourier transform of x
// y[j] = (1/n) sum x[k] e^(+2 pi i j k/n), so Inverse(Forward(x)) recovers x
//
// returns ErrInvalidLength if x is empty
//
// equivalent to scipy.fft.ifft(x) and MATLAB's ifft(x)
func Inverse(x []complex128) (y []complex128, err error) {
	if err = checkLength(len(x)); err != nil {
		return nil, err
	}
	// the inverse is the forward transform with conjugated twiddles
	// which is exatly conj(Forward(conj(x))): conjugation is exact in floating point
	y = conjugated(transform(conjugated(x)))
	scale := 1 / float64(len(x)) // multiplied, unlike pocketfft which divides
	for i := range y {
		y[i] *= complex(scale, 0)
	}
	return
}

// checkLength rejects lengths this package can't transform
// (for now, anything but power of two; mixed radix and Bluestein will change this when implemented)
func checkLength(n int) (err error) {
	if n == 0 {
		return fmt.Errorf("%w: cannot transform an empty signal", ErrInvalidLength)
	}
	if largestPrimeFactor(n) > maxRadix {
		return fmt.Errorf("%w: length %d has a prime factor above %d (not yet supported)", ErrInvalidLength, n, maxRadix)
	}
	return
}

// conjugated returns a new slice holding the complex conjuugate of each value in x
func conjugated(x []complex128) (c []complex128) {
	c = make([]complex128, len(x))
	for i, v := range x {
		c[i] = cmplx.Conj(v)
	}
	return
}

// transform picks the fastest algorithm for x's length
func transform(x []complex128) (y []complex128) {
	n := len(x)
	if isPowerOfTwo(n) {
		y = radix2(x)
	} else {
		y = mixedRadix(x, twiddleTable(n))
	}
	return
}
