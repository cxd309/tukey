package window

import (
	"fmt"
	"math"
)

// Tukey returns the n-point Tukey, or cosine-tapered, window:
// flat in the middle, with raised-cosine tapers at each end making up alpha of its length.
// alpha at or below 0 is Boxcar, and at or above 1 is Hann, infinities included.
//
//	w, err := window.Tukey(64, 0.25)
//
// Units:
//   - n is a number of samples
//   - alpha is the fraction of the window that tapers
//   - w is dimensionless
//
// Returns:
//   - ErrInvalidLength for a negative n.
//   - ErrInvalidParameter for a NaN alpha.
//   - ErrInvalidOption for an option other than Periodic.
//
// Equivalent to:
//   - scipy.signal.windows.tukey(n, alpha)
//   - MATLAB tukeywin(n, alpha), which has no periodic form
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

// tukey gives sample k of the symmetric Tukey window of length m, for 0 < alpha < 1:
// a raised-cosine taper over the first alpha·(m-1)/2 samples, then flat.
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
