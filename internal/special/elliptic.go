package special

import "math"

// EllipK is the complete elliptic integral of the first kind, K(m), for the parameter m = k² in [0, 1].
// K(1) is +Inf, and m outside [0, 1] gives NaN.
// For m near 1, where 1 - m has lost its digits, use EllipKM1.
func EllipK(m float64) (k float64) {
	return EllipKM1(1 - m)
}

// EllipKM1 is K(1 - p), taking the complementary parameter p directly,
// so a parameter near 1 keeps every digit of its distance from 1.
// EllipKM1(0) is +Inf, and p outside [0, 1] gives NaN.
// It's found by the arithmetic-geometric mean, same as scipy.special.ellipkm1(p):
//
//	K(1 - p) = π/(2·AGM(1, √p))
func EllipKM1(p float64) (k float64) {
	switch {
	case math.IsNaN(p) || p < 0 || p > 1:
		return math.NaN()
	case p == 0:
		return math.Inf(1)
	}
	a, b := 1.0, math.Sqrt(p)
	for range 64 { // it converges quadratically, in fewer than 10 for any p > 0
		if math.Abs(a-b) <= (2e-16)*a {
			break
		}
		a, b = (a+b)/2, math.Sqrt(a*b)
	}
	k = math.Pi / (a + b)
	return
}

// Jacobi returns the Jacobi elliptic functions sn, cn and dn of u for the parameter m,
// whose complement mc = 1 - m is given too, so neither loses digits when the other is near 0.
// m and mc must be in [0, 1] and sum to 1.
// It's by the descending Landen transformation,
// with the AGM's c_n found as c_n²/(4·a_(n+1)) rather than (a_n - b_n)/2, which cancels.
// dn is √(mc + m·cn²), a sum of positive terms, rather than Cephes's cn/cos(φ_1 - φ), which is 0/0 at u = K.
//
// For u in [0, K/2] each function is within a few ULP of its exact value while mc ≥ 1e-4,
// losing about two digits by mc = 1e-6, and four by 1e-9.
// In (K/2, K], cn is small, and its relative error grows as it shrinks,
// so a caller needing it there should use u′ = K - u, and JacobiInQuarter does:
//
//	cn(u) = √mc·sn(u′)/dn(u′)
//	dn(u) = √mc/dn(u′)
//	sn(u) = cn(u′)/dn(u′)
//
// It's scipy.special.ellipj(u, m)'s first three outputs, and MATLAB's ellipj(u, m),
// and more accurate than both as m approaches 1.
func Jacobi(u, m, mc float64) (sn, cn, dn float64) {
	switch {
	case mc == 0:
		sn, cn = math.Tanh(u), 1/math.Cosh(u)
		return sn, cn, cn
	case m == 0:
		sn, cn = math.Sincos(u)
		return sn, cn, 1
	}
	// a_n and c_n of the AGM of 1 and √mc, with c_0 = √m
	var a, c [16]float64
	a[0], c[0] = 1, math.Sqrt(m)
	b := math.Sqrt(mc)
	n := 0
	for n < len(a)-1 && c[n] > 1e-17*a[n] {
		a[n+1] = (a[n] + b) / 2
		c[n+1] = c[n] * c[n] / (4 * a[n+1])
		b = math.Sqrt(a[n] * b)
		n++
	}
	// the amplitude φ, from 2^n·a_n·u back down to n = 0
	phi := math.Ldexp(a[n]*u, n)
	for ; n > 0; n-- {
		phi = (math.Asin(c[n]*math.Sin(phi)/a[n]) + phi) / 2
	}
	sn, cn = math.Sincos(phi)
	dn = math.Sqrt(mc + m*cn*cn) // dn² = 1 - m·sn², as a sum of positive terms
	return
}

// CarlsonRF is Carlson's symmetric elliptic integral of the first kind,
// for x, y and z ≥ 0 with at most one of them 0:
//
//	RF(x, y, z) = ½·∫ dt/√((t + x)(t + y)(t + z)), for t from 0 to ∞
//
// It gives the incomplete integral of the first kind as:
//
//	F(φ | m) = sin φ·RF(cos²φ, 1 - m·sin²φ, 1)
//
// where the second argument can be formed without cancellation when m is near 1.
// It's by Carlson's duplication theorem, then his fifth-order series,
// which is within 1e-18 once the arguments are within 1e-3 of their mean,
// same as scipy.special.elliprf(x, y, z).
func CarlsonRF(x, y, z float64) (rf float64) {
	for range 64 { // each step shrinks the spread fourfold, so about 25 at most
		mean := (x + y + z) / 3
		if max(math.Abs(x-mean), math.Abs(y-mean), math.Abs(z-mean)) <= 1e-3*mean {
			break
		}
		sx, sy, sz := math.Sqrt(x), math.Sqrt(y), math.Sqrt(z)
		lambda := sx*sy + sy*sz + sz*sx
		x, y, z = (x+lambda)/4, (y+lambda)/4, (z+lambda)/4
	}
	mean := (x + y + z) / 3
	dx, dy := 1-x/mean, 1-y/mean
	dz := -(dx + dy)
	e2 := dx*dy - dz*dz
	e3 := dx * dy * dz
	rf = (1 - e2/10 + e3/14 + e2*e2/24 - 3*e2*e3/44) / math.Sqrt(mean)
	return
}

// ParameterFromNome is the parameter m whose nome, e^(-π·K(1 - m)/K(m)), is q, for q in [0, 1),
// by the theta series:
//
//	m = 16q·(Σ_(n≥0) q^(n(n+1))/(1 + 2·Σ_(n≥1) q^(n²)))^4
//
// The series needs at most 4 terms for q ≤ e^-π, where m ≤ ½.
// Above that, m is found as 1 - m′, where m′ is the parameter of the complementary nome q′,
// with ln q·ln q′ = π², and q′ < e^-π.
// m is then accurate, but 1 - m isn't, so a caller needing 1 - m should pass q′ itself.
// A q outside [0, 1) gives NaN.
// It's the same as mpmath.mfrom(q=q).
func ParameterFromNome(q float64) (m float64) {
	switch {
	case math.IsNaN(q) || q < 0 || q >= 1:
		return math.NaN()
	case q > nomeAtHalf: // m > ½, by the complementary nome, which is at most e^-π
		return 1 - parameterFromSmallNome(math.Exp(math.Pi*math.Pi/math.Log(q)))
	}
	m = parameterFromSmallNome(q)
	return
}

// IncompleteF is F(arctan(√tanSq) | 1 - m1),
// the incomplete elliptic integral of the first kind at a parameter near 1.
// It comes from Carlson's RF, whose second argument needs no subtraction:
//
//	1 - (1 - m1)·sin²φ = cos²φ + m1·sin²φ
func IncompleteF(tanSq, m1 float64) (f float64) {
	cosSq := 1 / (1 + tanSq)
	sinSq := tanSq / (1 + tanSq)
	f = math.Sqrt(sinSq) * CarlsonRF(cosSq, cosSq+m1*sinSq, 1)
	return
}

// JacobiInQuarter is Jacobi of u in [0, K], given uc = K - u too.
// Above K/2, cn is small and Jacobi loses its relative digits,
// so the functions come from uc instead, where Jacobi is accurate:
//
//	sn(u) = cn(uc)/dn(uc)
//	cn(u) = √mc·sn(uc)/dn(uc)
//	dn(u) = √mc/dn(uc)
func JacobiInQuarter(u, uc, m, mc float64) (sn, cn, dn float64) {
	if u <= uc {
		return Jacobi(u, m, mc)
	}
	snc, cnc, dnc := Jacobi(uc, m, mc)
	rootMC := math.Sqrt(mc)
	sn, cn, dn = cnc/dnc, rootMC*snc/dnc, rootMC/dnc
	return
}

// nomeAtHalf is e^-π, the nome where m = ½.
var nomeAtHalf = math.Exp(-math.Pi)

// parameterFromSmallNome is ParameterFromNome for q in [0, e^-π], where the series converges within 4 terms.
func parameterFromSmallNome(q float64) (m float64) {
	numerator, denominator := 1.0, 1.0 // the n = 0 terms
	for n := 1; n <= 8; n++ {
		term := math.Pow(q, float64(n*(n+1)))
		numerator += term
		denominator += 2 * math.Pow(q, float64(n*n))
		if term <= 1e-17*numerator {
			break
		}
	}
	ratio := numerator / denominator
	m = 16 * q * ratio * ratio * ratio * ratio
	return
}
