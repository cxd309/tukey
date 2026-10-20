package window

import "fmt"

// shape build an n-point window from value, which gives sample k of the
// symmetric window of length m, and handles what every window shares:
//   - n = 0 is empty and n=1 is [1], before any formula, same as SciPy
//   - periodic is the symmetric window of length n+1 without its last point
//   - symmetric is exact: value is evaluated for the first half only
//     including the middle sample when m is odd, and mirrored
//
// returns ErrInvalidLength if n < 0
func shape(n int, periodic bool, value func(k, m int) float64) (w []float64, err error) {
	switch {
	case n < 0:
		return nil, fmt.Errorf("%w: n must be >= 0, got %d", ErrInvalidLength, n)
	case n == 0:
		return []float64{}, nil
	case n == 1:
		return []float64{1}, nil
	}
	m := n
	if periodic {
		m = n + 1
	}
	w = make([]float64, n)
	for k := range (m + 1) / 2 {
		v := value(k, m)
		w[k] = v
		if mirror := m - 1 - k; mirror < n {
			w[mirror] = v
		}
	}
	return
}
