package window

import (
	"fmt"
	"math"
)

// Taylor returns the n-point Taylor window, for radar and array processing:
// its first nbar-1 sidelobes are nearly level at sll dB below the main lobe,
// and the rest fall away.
// Normalised to 1 at the centre unless Unnormalised.
// nbar = 1 is Boxcar
//
// returns ErrInvalidLength if n < 0
// returns ErrInvalidParameter if nbar < 1 or unless sll > 0
//
// equivalent to scipy.signal.windows.taylor(n, nbar, sll)
// equivalent to MATLAB taylorwin(n, nbar, -sll) is Taylor with Unnormalised, for n > 1
func Taylor(n, nbar int, sll float64, opts ...Option) (w []float64, err error) {
	cfg, err := newConfig("Taylor", opts, optionSet{optionPeriodic: true, optionUnnormalised: true})
	if err != nil {
		return nil, err
	}
	if nbar < 1 {
		return nil, fmt.Errorf("%w: Taylor's nbar must be >= 1, got %d", ErrInvalidParameter, nbar)
	}
	if err = positive("Taylor", "sll", sll); err != nil {
		return nil, err
	}
	f := taylorCoefficients(nbar, sll)
	value := func(k, m int) float64 { return taylorSum(f, float64(k), m) }
	if !cfg.unnormalised {
		value = func(k, m int) float64 {
			return taylorSum(f, float64(k), m) / taylorSum(f, float64(m-1)/2, m)
		}
	}
	w, err = shape(n, cfg.periodic, value)
	return
}

// taylorCoefficients is F[1..nbar-1], the weights of Taylor's cosine terms,
// which depend only on nbar and the sidelobe level.
// f[0] is unused, so f[m] is F_m
func taylorCoefficients(nbar int, sll float64) (f []float64) {
	a := math.Acosh(math.Pow(10, sll/20)) / math.Pi
	s2 := float64(nbar*nbar) / (a*a + (float64(nbar)-0.5)*(float64(nbar)-0.5))
	f = make([]float64, nbar)
	for m := 1; m < nbar; m++ {
		mm := float64(m * m)
		numer, denom := 1.0, 2.0
		for i := 1; i < nbar; i++ {
			numer *= 1 - mm/s2/(a*a+(float64(i)-0.5)*(float64(i)-0.5))
			if i != m {
				denom *= 1 - mm/float64(i*i)
			}
		}
		if m%2 == 0 {
			numer = -numer
		}
		f[m] = numer / denom
	}
	return
}

// taylorSum is the Taylor window of length m at position x:
// 1 plus twice the sum of F_m cos(2 pi m (x - m/2 + 1/2) / m).
// at x = (m-1)/2 it's the centre
func taylorSum(f []float64, x float64, m int) (v float64) {
	v = 1
	for j := 1; j < len(f); j++ {
		v += 2 * f[j] * math.Cos(2*math.Pi*float64(j)*(x-float64(m)/2+0.5)/float64(m))
	}
	return
}
