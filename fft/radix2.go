package fft

import "math/bits"

// radix 2 is the iterative Cooley-Tukey FFT for a power of two length
// the input is put in bit-reversed order
// then log2(n) stages each combine pairs of transforms of size half
// into transforms of size 2*half with n/2 butterflies
//
// X[k] 		= E[k] + w^k·O[k]
// X[k+half] 	= E[k] - w^k·O[k]
func radix2(x []complex128) (y []complex128) {
	n := len(x)
	y = bitReversed(x)

	// every stage's twiddles are a subset of e^(-2 pi i j/n) for j < n/2
	table := make([]complex128, n/2)
	for j := range table {
		table[j] = twiddle(j, n)
	}

	for size := 2; size <= n; size *= 2 {
		half, stride := size/2, n/size
		for start := 0; start < n; start += size {
			for j := range half {
				u := y[start+j]
				t := table[j*stride] * y[start+j+half]
				y[start+j] = u + t
				y[start+j+half] = u - t
			}
		}
	}
	return
}

// bitReversed returns a copy of x with element i moved to the index whose binary digits are i's reversed
// the order the butterfly stages need
func bitReversed(x []complex128) (y []complex128) {
	logN := bits.Len(uint(len(x))) - 1
	y = make([]complex128, len(x))
	for i, v := range x {
		y[bits.Reverse(uint(i))>>(bits.UintSize-logN)] = v
	}
	return
}

// isPowerOfTwo reports whether n is a positive power of two
func isPowerOfTwo(n int) (ok bool) {
	ok = n > 0 && n&(n-1) == 0
	return
}
