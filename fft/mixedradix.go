package fft

// maxRadix is the largest prime factor whose butterfly is a direct sum
// larger primes use Bluestein
const maxRadix = 13

// mixedRadix is the recursive Cooley-Tukey FFT for lengths whose prime factors are all at most maxRadix
// with n = pm and p the smallest prime factor,
// x splits into p interleaved subsequences x[r], x[r+p], x[r+2p], ... of length m
// each transformed recursively as Y_r, output k=q+sm then combines them
//
//	X[q + sm] = sum_r (w_n^(rq) Y_r[q]) w_p^(rs)
//
// a twiddle between levels, then a p-point DFT (the butterfly)
// table holds e^(-2 pi i j/N) for the top-level length N,
// so w_n^j = table[j N/n]
func mixedRadix(x []complex128, table []complex128) (y []complex128) {
	n := len(x)
	if n == 1 {
		y = []complex128{x[0]}
		return
	}
	p := smallestPrimeFactor(n)
	m := n / p
	stride := len(table) / n // table[j*stride] is w_n^j

	// n is prime: there's nothing to split into
	// so the whole tranform is one p-point butterfly (each subsequence would be a single sample)
	if m == 1 {
		y = butterfly(x, table, stride)
		return
	}

	// transform each interleaved subsequence
	transformed := make([][]complex128, p)
	sub := make([]complex128, m)
	for r := range p {
		for j := range m {
			sub[j] = x[r+j*p]
		}
		transformed[r] = mixedRadix(sub, table)
	}

	//combine: for each q, twiddle the p values, then a p-point DFT across time
	y = make([]complex128, n)
	terms := make([]complex128, p)
	for q := range m {
		for r := range p {
			terms[r] = table[r*q*stride] * transformed[r][q]
		}
		out := butterfly(terms, table, m*stride)
		for s := range p {
			y[q+s*m] = out[s]
		}
	}
	return
}

// smallestPrimeFactor returns n's smallest prime factor or n itself if n is prime (n<=2)
func smallestPrimeFactor(n int) (p int) {
	for p = 2; p*p <= n; p++ {
		if n%p == 0 {
			return
		}
	}
	return n
}

// largestPrimeFactor returns n's largest prime factor, or 1 for n=1
func largestPrimeFactor(n int) (p int) {
	p = 1
	for f := 2; f*f <= n; f++ {
		for n%f == 0 {
			p, n = f, n/f
		}
	}
	if n > 1 {
		p = n
	}
	return
}

// butterfly is the p-point DFT at the heart of each mixed-radix combination:
// a direct sum for primes up to maxRadix,
// Bluestein's algorithm for larger ones, so every length stays O(n log n).
// w_p^j is table[j*wStep].
func butterfly(terms, table []complex128, wStep int) (out []complex128) {
	p := len(terms)
	if p > maxRadix {
		out = bluestein(terms)
		return
	}
	out = make([]complex128, p)
	for s := range p {
		var sum complex128
		for r := range p {
			sum += terms[r] * table[(r*s%p)*wStep]
		}
		out[s] = sum
	}
	return
}
