package dsp

import (
	"math"

	"github.com/cxd309/tukey/internal/poly"
	"github.com/cxd309/tukey/internal/special"
)

// ellipticPrototype returns the analog elliptic lowpass prototype:
// equal ripple of rippleDB up to 1 rad/s,
// and equal ripple at least attenuationDB down from the stopband edge 1/√m rad/s,
// with zeros on the imaginary axis.
// An even order's gain is set so the passband's peak is 0dB, which puts DC at -rippleDB, as for Cheby1.
//
// The design is Orfanidis's, which scipy ellipap also follows,
// but every elliptic parameter is carried with its complement,
// so a sharp design whose m is near 1 keeps the digits of 1 - m:
//   - the degree equation N·K`(m)/K(m) = K`(m1)/K(m1), solved by both nomes:
//     ln q = ln q1 / N for m, and ln q` = N·ln q1` for 1 − m
//   - the inverse sc, sc(r | 1 − m1) = 1/ε, by Carlson's RF, with K`(m1) − r
//     found the same way rather than by subtracting
//   - sn, cn and dn above half a quarter period from their values at the
//     complementary argument, where Jacobi is accurate
//
// Equivalent to scipy.signal.ellipap.
func ellipticPrototype(order int, rippleDB, attenuationDB float64) (proto analogZPK) {
	epsilonSq := math.Expm1(rippleDB * math.Ln10 / 10)
	if order == 1 {
		pole := -1 / math.Sqrt(epsilonSq)
		proto.Poles = []complex128{complex(pole, 0)}
		proto.Gain = -pole
		return
	}
	stopSq := math.Expm1(attenuationDB * math.Ln10 / 10)
	m1 := epsilonSq / stopSq // the specification's modulus, squared
	k1, k1Prime := special.EllipKM1(1-m1), special.EllipKM1(m1)

	// the degree equation, by the nomes q = q1^(1/N) and q` = q1`^N
	n := float64(order)
	m := special.ParameterFromNome(math.Exp(-math.Pi * k1Prime / (n * k1)))
	mc := special.ParameterFromNome(math.Exp(-math.Pi * n * k1 / k1Prime))
	k, kPrime := special.EllipKM1(mc), special.EllipKM1(m)

	// v0 = K`·r/K`(m1), where sc(r | 1-m1)=1/ε, so tan phi = 1/ε.
	// its complement K`- v0 is K`·r̄K`(m1), where tan psi = 1/(√m1·tan phi) = √stopSq
	r := incompleteF(1/epsilonSq, m1)
	rBar := incompleteF(stopSq, m1)
	sv, cv, dv := jacobiInQuarter(kPrime*r/k1Prime, kPrime*rBar/k1Prime, mc, m)

	// the j-th of N points across the quarter period K: j=0,
	// in an odd order only, is the real pole and a zero at infinity
	for j := 1 - order%2; j < order; j += 2 {
		s, c, d := jacobiInQuarter(float64(j)*k/n, float64(order-j)*k/n, m, mc)
		// 1 − (d·sv)², as a sum of positive terms
		denominator := cv*cv + m*(s*sv)*(s*sv)
		pole := -complex(c*d*sv*cv, s*dv) / complex(denominator, 0)
		if j == 0 {
			proto.Poles = append(proto.Poles, complex(real(pole), 0))
			continue
		}
		zero := complex(0, 1/(math.Sqrt(m)*s))
		proto.Zeros = append(proto.Zeros, zero, -zero)
		proto.Poles = append(proto.Poles, pole, complex(real(pole), -imag(pole)))
	}
	proto.Gain = real(poly.Prod(poly.MapRoots(proto.Poles, negate)) / poly.Prod(poly.MapRoots(proto.Zeros, negate)))
	if order%2 == 0 {
		proto.Gain *= math.Pow(10, -rippleDB/20)
	}
	return
}

// incompleteF is F(arctan(√tanSq) | 1-m1),
// the incomplete elliptic integral of the first kind at a parameter near 1,
// from Carlson's RF, whose second argment
// 1- (1- m1)·sin^2 phi = cos^2 phi + m1·sin^2 phi needs no subtraction
func incompleteF(tanSq, m1 float64) (f float64) {
	cosSq := 1 / (1 + tanSq)
	sinSq := tanSq / (1 + tanSq)
	f = math.Sqrt(sinSq) * special.CarlsonRF(cosSq, cosSq+m1*sinSq, 1)
	return
}

// jacobiInQuarter is special.Jacobi of u in [0,K], given uc=K- u too.
// Above K/2, cn is small and Jacobi loses its relative digits,
// so the functions come from uc instead, where Jacobi is accurate:
//   - sn(u)= cn(uc)/ dn(uc),
//   - cn(u)= √mc·sn(uc)/ dn(uc), and,
//   - dn(u)= √mc/ dn(uc).
func jacobiInQuarter(u, uc, m, mc float64) (sn, cn, dn float64) {
	if u <= uc {
		return special.Jacobi(u, m, mc)
	}
	snc, cnc, dnc := special.Jacobi(uc, m, mc)
	rootMC := math.Sqrt(mc)
	sn, cn, dn = cnc/dnc, rootMC*snc/dnc, rootMC/dnc
	return
}

// Ellip designs a digital elliptic (Cauer) filter for band and returns ZPK.
// It has equal ripple of rippleDB in the passband and equal ripple at least attenuationDB down in the stopband,
// the sharpest transition of any filter of its order.
// A band edge is where the gain first drops below -rippleDB, as for Cheby1.
//
//	f, err:= dsp.Ellip(4, 1, 40, dsp.Lowpass(0.3))
//
// Returns:
//   - ErrInvalidOrder if order < 1
//   - ErrInvalidBand for a band Butter rejects
//   - ErrInvalidParameter if rippleDB or attenuationDB isn't positive and finite,
//     or attenuationDB isn't greater than rippleDB
//
// Equivalent to:
//   - scipy.signal.ellip(order, rippleDB, attenuationDB, wn, btype, output='zpk')
//   - matlab ellip(order, rippleDB, attenuationDB, wn, ftype)
func Ellip(order int, rippleDB, attenuationDB float64, band Band) (f ZPK, err error) {
	if err = checkRippleAttenuation(rippleDB, attenuationDB); err != nil {
		return f, err
	}
	return designDigital(order, band, func() analogZPK { return ellipticPrototype(order, rippleDB, attenuationDB) })
}
