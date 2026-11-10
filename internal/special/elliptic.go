package special

import "math"

// EllipK is the complete elliptic integral of the first kind, K(m),
// for the parameter m=k^2 in [0,1].
// K(1) is +Inf, and m outside [0,1] gives NaN.
// For m near 1, where 1-m has lost its digits, use EllipKM1.
func EllipK(m float64) (k float64) {
	return EllipKM1(1 - m)
}

// EllipKM1 is K(1-p), taking the complementary parameter p directly,
// so a parameter near 1 keeps every digit of its distance from 1.
// EllipKM1(0) is +Inf, and p outside [0,1] gives NaN.
//
// K(1-p)=pi / (2·AGM(1, sqrt(p))), by the arithmetic-geometric mean.
//
// Equivalent to scipy.special.ellipkm1(p)
func EllipKM1(p float64) (k float64) {
	switch {
	case math.IsNaN(p) || p < 0 || p > 1:
		return math.NaN()
	case p == 0:
		return math.Inf(1)
	}
	a, b := 1.0, math.Sqrt(p)
	for range 64 { // quadratic convergence: fewer than 10 for any p > 0
		if math.Abs(a-b) <= (2e-16)*a {
			break
		}
		a, b = (a+b)/2, math.Sqrt(a*b)
	}
	k = math.Pi / (a + b)
	return
}

// Jacobi returns the Jacobi elliptic functions sn, cn, and dn of use for the parameter m,
// whose complement mc = 1-m is given too, so neither loses digits when the other is near 0.
// m and mc must be in [0, 1] and sum to 1.
//
// By the descending Landen transformation, the AGM's c_n found as
// c_n^2 / (4·a_{n+1}) rather than (a_n - b_n)/2, which cancels.
// dn is sqrt(mc + m·cn^2), a sum of positive terms,
// rather than Cephes's cn/cos(phi_1 - phi), which is 0/0 at u=K.
//
// For u in [0, K/2] each function is within a few ULP of its exact value while mc >= 1e-4,
// losing about two digitis by mc = 1e-6 and found by 1e-9.
// In (K/2, K], cn is small, and its relative error grows as it shrinks:
// a caller needing it there should use u′ = K-u, where
// 	- cn(u) = √mc·sn(u′)/dn(u′),
// 	- dn(u) = √mc/dn(u′), and,
// 	- sn(u) = cn(u′)/dn(u′).

// Equivalent to:
//   - scipy.special.ellipj(u, m) (first 3 outputs)
//   - matlab ellipj(u, m)
//
// and is more accurate than both as m approaches 1
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
	// the amplitude phi, from 2^n·a_n·u back down to n = 0
	phi := math.Ldexp(a[n]*u, n)
	for ; n > 0; n-- {
		phi = (math.Asin(c[n]*math.Sin(phi)/a[n]) + phi) / 2
	}
	sn, cn = math.Sincos(phi)
	dn = math.Sqrt(mc + m*cn*cn) // dn^2 = 1-m·sn^2, as a sum of positive terms
	return
}

// CarlsonRF is Carlson's symmetric elliptic integral of the first kind,
// RF(x,y,z) = 1/2 integral e^inf dt / sqrt((t+x)(t+z)), for x,y,z >=0 with at most one of them 0.
// It gives the incomplete integral of the first kind as
// F(phi | m) = sin phi · RF(cos^2 phi, 1-m·sin^2 phi, 1),
// where the second argument can be formed without cancellation when m is near 1.
//
// By Carlson's duplication theorum, then his fifth-order series,
// which is within 1e-18 once the arguments are within 1e-3 of their mean.
//
// Equivalent to:
//   - scipy.special.elliprf(x,y,z)
func CarlsonRF(x, y, z float64) (rf float64) {
	for range 64 { // each step shrinks the spread fourfold: about 25 at most
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

// nomeAtHalf is e^-pi, the nome where m=1/2
var nomeAtHalf = math.Exp(-math.Pi)

// ParameterFromNome is the parameter m whose nome q = exp(-pi·K(1-m)/K(m))
// is q, for q in [0,1), by the theta series
// m = 16q·(Σ_{n≥0} q^{n(n+1)} / (1 + 2·Σ_{n≥1} q^{n²}))⁴.
// The series needs at most 4 terms for q ≤ e^-pi, where m ≤ 1/2. Above
// that, m is found as 1 - m′, where m′ is the parameter of the
// complementary nome q′, with ln q · ln q′ = pi^2, and q′ < e^-pi.
// m is then accurate, but 1 - m is not: a caller needing 1 - m should
// pass q′ itself.
// q outside [0,1) gives NaN.
//
// Equivalent to:
//   - mpmath.mfrom(q=q).
func ParameterFromNome(q float64) (m float64) {
	switch {
	case math.IsNaN(q) || q < 0 || q >= 1:
		return math.NaN()
	case q > nomeAtHalf: // m > ½: by the complementary nome, at most e^−π
		return 1 - parameterFromSmallNome(math.Exp(math.Pi*math.Pi/math.Log(q)))
	}
	m = parameterFromSmallNome(q)
	return
}

// parameterFromSmallNome is ParameterFromNome for q in [0, e^−π], where the
// series converges within 4 terms.
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
