package dsp

import (
	"fmt"
	"math"
	"math/cmplx"

	"github.com/cxd309/tukey/internal/poly"
)

// rippleFactor is ε= sqrt(10^(dB/10)- 1), by Expm1 so a small dB keeps its digits
func rippleFactor(dB float64) (epsilon float64) {
	epsilon = math.Sqrt(math.Expm1(dB * math.Ln10 / 10))
	return
}

// checkDB returns ErrInvalidParameter unless dB is positive and finite
func checkDB(name string, dB float64) (err error) {
	if !(dB > 0) || math.IsInf(dB, 1) {
		err = fmt.Errorf("%w: %s must be positive and finite, got %g", ErrInvalidParameter, name, dB)
	}
	return
}

// chebyshev1Prototype returns the analog Chebyshev type I lowpass prototype:
// equal ripple of rippleDB up to 1 rad/s, poles on an ellipse and no finite zeros.
// An even order's gain is set so the passband's peak is 0dB, which puts DC at -rippleDB.
//
// Equivalent to scipy.signal.cheb1ap
func chebyshev1Prototype(order int, rippleDB float64) (proto analogZPK) {
	mu := math.Asinh(1/rippleFactor(rippleDB)) / float64(order)
	proto.Poles = make([]complex128, order)
	for k := range order {
		proto.Poles[k] = -cmplx.Sinh(complex(mu, prototypeAngle(order, k)))
	}
	proto.Gain = real(poly.Prod(poly.MapRoots(proto.Poles, negate)))
	if order%2 == 0 {
		proto.Gain *= math.Pow(10, -rippleDB/20)
	}
	return
}

// Cheby1 designs a digital Chebyshev type I filter for band and returns ZPK.
// It has equal ripple of rippleDB in the passband and falls monotonically after it,
// a sharper transition than Butter's for the same order.
// A band edge is where the gain first drops below -rippleDB, not -3 dB
//
//	f, err := dsp.Cheby1(4, 1, dsp.Lowpass(0.3))
//
// Returns:
//   - ErrInvalidOrder if order < 1
//   - ErrInvalidBand for a band Butter rejects
//   - ErrInvalidParameter if rippleDB isn't positive and finite
//
// Equivalent to:
//   - scipy.signal.cheby1(order, rippleDB, wn, btype, output='zpk')
//   - matlab cheby1(order, rippleDB, wn, ftype)
func Cheby1(order int, rippleDB float64, band Band) (f ZPK, err error) {
	if err = checkDB("rippleDB", rippleDB); err != nil {
		return f, err
	}
	return designDigital(order, band, func() analogZPK { return chebyshev1Prototype(order, rippleDB) })
}

// chebyshev2Prototype returns the analog Chebyshev type II lowpass prototype:
// flat in the passband, with equal ripple at least attenuationDB down from 1 rad/s.
// Zeros are on the imaginary axis, and an odd order's middle zero is at infinity.
//
// Equivalent to scipy.signal.cheb2ap.
func chebyshev2Prototype(order int, attenuationDB float64) (proto analogZPK) {
	mu := math.Asinh(rippleFactor(attenuationDB)) / float64(order)
	proto.Poles = make([]complex128, order)
	for k := range order {
		proto.Poles[k] = -1 / cmplx.Sinh(complex(mu, prototypeAngle(order, k)))
		if m := -order + 1 + 2*k; m != 0 {
			proto.Zeros = append(proto.Zeros, complex(0, 1/math.Sin(prototypeAngle(order, k))))
		}
	}
	proto.Gain = real(poly.Prod(poly.MapRoots(proto.Poles, negate)) / poly.Prod(poly.MapRoots(proto.Zeros, negate)))
	return
}

// Cheby2 designs a digital Chebyshev type II filter for band and returns ZPK.
// It's flat in the passband, with equal ripple in the stopband at least attenuationDB down.
// A band edge is the start of the stopband, where the gain first reaches -attenuationDB,
// not the end of the passband as for Butter and Cheby1.
//
//	f, err := dsp.Cheby2(4, 40, dsp.Lowpass(0.3))
//
// Returns:
//   - ErrInvalidOrder if order < 1
//   - ErrInvalidBand for a band Butter rejects
//   - ErrInvalidParameter if attenuationDB isn't positive and finite
//
// Equivalent to:
//   - scipy.signal.cheby2(order, attenuationDB, wn, btype, output='zpk')
//   - matlab cheby2(order, attenuationDB, wn, ftype)
func Cheby2(order int, attenuationDB float64, band Band) (f ZPK, err error) {
	if err = checkDB("attenuationDB", attenuationDB); err != nil {
		return f, err
	}
	return designDigital(order, band, func() analogZPK { return chebyshev2Prototype(order, attenuationDB) })
}
