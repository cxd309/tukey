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
func mixedRadix(x []complex128) (y []complex128) {
	n := len(x)
	// a large prime is a single Bluestein butterfly
	// doesn't need the twiddle table or the scratch buffer so don't build
	if n > maxRadix && smallestPrimeFactor(n) == n {
		y = bluestein(x)
		return
	}
	y = make([]complex128, n)
	plan := mixedRadixPlan{
		table: twiddleTable(n),
		terms: make([]complex128, max(largestPrimeFactor(n), 4)),
	}
	plan.transform(x, 1, y)
	return
}

// mixedRadixPlan is what every level of one transform shares
type mixedRadixPlan struct {
	table []complex128 // e^(-2 pi i·j/N) for the top-level length N
	terms []complex128 // scratch for one butterfly's inputs, sized for the largest factor
}

// transform writes the DFT of in[0], in[step], in[2·step], ... (len(out) values) to out
//
// it works in place in out:
// sub-transform r goes to block out[rm : (r+1)m],
// each combination reads the p values at q, q+m, ..., q+(p-1)m,
// which are exactly the positions its outputs go to
// subsequence r is read straight from in, starting at in[r·step] with step·p between values, so nothing is copied
func (plan mixedRadixPlan) transform(in []complex128, step int, out []complex128) {
	n := len(out)
	if n == 1 {
		out[0] = in[0]
		return
	}
	p := nextFactor(n)
	m := n / p
	stride := len(plan.table) / n // plan.table[j*stride] is w_n^j
	terms := plan.terms[:p]

	// n is a single butterfly: there's nothing to split into, so the whole transform is one butterfly
	if m == 1 {
		for r := range p {
			terms[r] = in[r*step]
		}
		butterfly(terms, plan.table, stride, out, 1)
		return
	}

	for r := range p {
		plan.transform(in[r*step:], step*p, out[r*m:(r+1)*m])
	}
	for q := range m {
		for r := range p {
			terms[r] = plan.table[r*q*stride] * out[r*m+q]
		}
		butterfly(terms, plan.table, m*stride, out[q:], m)
	}
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

// butterfly writes the p-point DFT of terms to dst[0], dst[step], ..., dst[(p-1)*step]:
// dedicated versions for 2, 3, 4 and 5, a direct sum for other primes up to
// maxRadix, and Bluestein's algorithm for larger ones. w_p^j is table[j*wStep]
// it only reads terms, so dst can hold the values terms was built from.
func butterfly(terms, table []complex128, wStep int, dst []complex128, step int) {
	switch p := len(terms); {
	case p == 2:
		butterfly2(terms, dst, step)
	case p == 3:
		butterfly3(terms, dst, step)
	case p == 4:
		butterfly4(terms, dst, step)
	case p == 5:
		butterfly5(terms, dst, step)
	case p > maxRadix:
		for s, v := range bluestein(terms) {
			dst[s*step] = v
		}
	default:
		directButterfly(terms, table, wStep, dst, step)
	}
}

// directButterfly is the p-point DFT as a direct sum, for primes with no dedicated butterfly
func directButterfly(terms, table []complex128, wStep int, dst []complex128, step int) {
	p := len(terms)
	for s := range p {
		var sum complex128
		for r := range p {
			sum += terms[r] * table[(r*s%p)*wStep]
		}
		dst[s*step] = sum
	}
}

// nextFactor is the radix n is split by: 4 when it divides n,
// since one radix-4 level costs less than two radix-2 ones,
// otherwise n's smallest prime factor
func nextFactor(n int) (p int) {
	if n%4 == 0 {
		return 4
	}
	return smallestPrimeFactor(n)
}
