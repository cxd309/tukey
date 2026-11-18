package fft

// The constants for the 3- and 5-point butterflies, as high-precision literals.
// The compiler rounds each to the nearest float64,
// slightly more accurate than evaluating sin and cos at run time.
const (
	sin60  = 0.86602540378443864676372317075293618347140262690519 // sin(π/3) = √3/2
	cos72  = 0.30901699437494742410229341718281905886436644335    // cos(2π/5)
	cos144 = -0.80901699437494742410229341718281905886436644335   // cos(4π/5)
	sin72  = 0.95105651629515357211643933337938214340569863412575 // sin(2π/5)
	sin144 = 0.58778525229247312916870595463907276859765243764314 // sin(4π/5)
)

// butterfly2 is the 2-point DFT:
//
//	X0 = t0 + t1, X1 = t0 - t1
func butterfly2(t, dst []complex128, step int) {
	dst[0], dst[step] = t[0]+t[1], t[0]-t[1]
}

// butterfly3 is the 3-point DFT, with w = e^(-2πi/3) = -1/2 - i·√3/2:
//
//	X0 = t0 + (t1 + t2)
//	X1 = t0 - (t1 + t2)/2 - i·(√3/2)·(t1 - t2)
//	X2 = t0 - (t1 + t2)/2 + i·(√3/2)·(t1 - t2)
func butterfly3(t, dst []complex128, step int) {
	sum, diff := t[1]+t[2], t[1]-t[2]
	mid := t[0] - scaledBy(0.5, sum)
	rot := minusI(scaledBy(sin60, diff))
	dst[0] = t[0] + sum
	dst[step], dst[2*step] = mid+rot, mid-rot
}

// butterfly4 is the 4-point DFT, with w = e^(-2πi/4) = -i, so every rotation only swaps and negates:
//
//	X0 = (t0 + t2) + (t1 + t3)    X2 = (t0 + t2) - (t1 + t3)
//	X1 = (t0 - t2) - i·(t1 - t3)  X3 = (t0 - t2) + i·(t1 - t3)
func butterfly4(t, dst []complex128, step int) {
	a, b := t[0]+t[2], t[0]-t[2]
	c, d := t[1]+t[3], t[1]-t[3]
	rot := minusI(d)
	dst[0], dst[2*step] = a+c, a-c
	dst[step], dst[3*step] = b+rot, b-rot
}

// butterfly5 is the 5-point DFT.
// It pairs the inputs symmetrically about t0,
// so each pair adds a cosine to the real-like part and a sine to the rotated part:
//
//	a1 = t1 + t4, b1 = t1 - t4, a2 = t2 + t3, b2 = t2 - t3
//	X0 = t0 + a1 + a2
//	X1, X4 = t0 + cos72·a1 + cos144·a2 ∓ i·(sin72·b1 + sin144·b2)
//	X2, X3 = t0 + cos144·a1 + cos72·a2 ∓ i·(sin144·b1 - sin72·b2)
func butterfly5(t, dst []complex128, step int) {
	a1, b1 := t[1]+t[4], t[1]-t[4]
	a2, b2 := t[2]+t[3], t[2]-t[3]

	r1 := t[0] + scaledBy(cos72, a1) + scaledBy(cos144, a2)
	r2 := t[0] + scaledBy(cos144, a1) + scaledBy(cos72, a2)
	i1 := minusI(scaledBy(sin72, b1) + scaledBy(sin144, b2))
	i2 := minusI(scaledBy(sin144, b1) - scaledBy(sin72, b2))

	dst[0] = t[0] + a1 + a2
	dst[step], dst[4*step] = r1+i1, r1-i1
	dst[2*step], dst[3*step] = r2+i2, r2-i2
}

// scaledBy multiplies z by the real number k:
// two real multiplications, where a complex multiplication by k + 0i would do four.
func scaledBy(k float64, z complex128) (scaled complex128) {
	scaled = complex(k*real(z), k*imag(z))
	return
}

// minusI multiplies z by -i exactly, a swap and a negation:
//
//	-i·(a + bi) = b - ai
func minusI(z complex128) (rotated complex128) {
	rotated = complex(imag(z), -real(z))
	return
}
