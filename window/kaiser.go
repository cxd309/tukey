package window

import (
	"fmt"
	"math"

	"github.com/cxd309/tukey/internal/special"
)

// Kaiser returns the n-point Kaiser windows with shape beta:
// I0(beta sqrt(1 - x^2)) / I0(beta), for x from -1 to 1
// beta 0 is Boxcar, and larger beta traders a wider main love for lower
// sidelobes: about 5 is like Hamming, 8.6 like Blackman, 14 for very low sidelobes
// a negative beta is the same window as |beta|, as I0 is even
//
// returns ErrInvalidLength if n<0, and ErrInvalidParameter if beta is NaN or +=Inf
//
// equivalent to scipy.signal.windows.kaiser(n, beta),
// equivalent to MATLAB kaiser(n, beta)
func Kaiser(n int, beta float64, opts ...Option) (w []float64, err error) {
	if math.IsNaN(beta) || math.IsInf(beta, 0) {
		return nil, fmt.Errorf("%w: Kaiser's beta must be finite, got %g", ErrInvalidParameter, beta)
	}
	peak := special.I0(beta)
	w, err = newWindow("Kaiser", n, opts, func(k, m int) float64 {
		half := float64(m-1) / 2
		x := (float64(k) - half) / half
		return special.I0(beta*math.Sqrt(1-x*x)) / peak
	})
	return
}

// KaiserBesselDerived returns the n-point Kaiser-Bessel-derived (KBD) window,
// the MDCT window of audio codecs such as AAC and Vorbis:
// the square root of the running sum of a Kaiser window,
// so that w[k]^2 + w[k+n/2]^2 = 1
// (the Princeton-Bradley condition, for perfect reconstruction)
// it's defined only for even n, and only symmetric
//
// returns ErrInvalidLength if n < 0 or n is odd (including n = 1),
// returns ErrInvalidParameter if beta is NaN or +=Inf, or given any option
//
// equivalent to scipy.signal.windows.kaiser_bessel_derived(n, beta),
// no MATLAB equivalent
func KaiserBesselDerived(n int, beta float64, opts ...Option) (w []float64, err error) {
	if _, err = newConfig("KaiserBesselDerived", opts, optionSet{}); err != nil {
		return nil, err
	}
	switch {
	case n < 0:
		return nil, fmt.Errorf("%w: n must be >=0, got %d", ErrInvalidLength, n)
	case n%2 == 1:
		return nil, fmt.Errorf("%w: KaiserBesselDerived needs an even n, got %d", ErrInvalidLength, n)
	case n == 0:
		return []float64{}, nil
	}
	kaiser, err := Kaiser(n/2+1, beta)
	if err != nil {
		return nil, err
	}
	total := 0.0
	for _, v := range kaiser {
		total += v
	}
	w = make([]float64, n)
	running := 0.0
	for k := range n / 2 {
		running += kaiser[k]
		v := math.Sqrt(running / total)
		w[k], w[n-1-k] = v, v
	}
	return
}
