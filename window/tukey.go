package window

import (
	"fmt"
	"math"
)

// Tukey returns the n-point tukey (cosine-tapered) window:
// flat in the middle with raised-cosine tapers at each end
// making up alpha of its length
//
// alpha 0 or less is Boxcar, and 1 or more is Hann
//
// returns ErrInvalidLength if n < 0, and ErrInvalidParameter if alpha is NaN
//
// equivalent to scipy.signal.windows.tukey(n, alpha) and MATLAB's tukeywin(n, alpha)
func Tukey(n int, alpha float64, opts ...Option) (w []float64, err error) {
	switch {
	case math.IsNaN(alpha):
		return nil, fmt.Errorf("%w: Tukey's alpha is NaN", ErrInvalidParameter)
	case alpha <= 0:
		w, err = newWindow("Tukey", n, opts, boxcar)
	case alpha >= 1:
		w, err = newWindow("Tukey", n, opts, cosineSum(generalHamming(0.5)))
	default:
		w, err = newWindow("Tukey", n, opts, tukey(alpha))
	}
	return
}

// tukey gives sample k of the symmetric Tukey window of length m, for 0<alpha<1:
// a raised-cosine taper over the first alpha (m-1)/2 samples, then flat
func tukey(alpha float64) (value func(k, m int) float64) {
	value = func(k, m int) (v float64) {
		width := int(math.Floor(alpha * float64(m-1) / 2))
		if k > width {
			return 1
		}
		v = 0.5 * (1 + math.Cos(math.Pi*(-1+2*float64(k)/(alpha*float64(m-1)))))
		return
	}
	return
}
