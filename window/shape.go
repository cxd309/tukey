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
	if w, short, err := shortWindow(n); short {
		return w, err
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

// sampled builds an n-point window that isn't symmetric, from value,
// which gives sample i directly, lengths are same as shape
func sampled(n int, value func(i int) float64) (w []float64, err error) {
	if w, short, err := shortWindow(n); short {
		return w, err
	}
	w = make([]float64, n)
	for i := range w {
		w[i] = value(i)
	}
	return
}

// newWindow builds the n-point window name from value (see shape)
// accepting only Periodic
func newWindow(name string, n int, opts []Option, value func(k, m int) float64) (w []float64, err error) {
	cfg, err := newConfig(name, opts, optionSet{optionPeriodic: true})
	if err != nil {
		return nil, err
	}
	w, err = shape(n, cfg.periodic, value)
	return
}

// shortWindow handles the lengths every window shares, before any formula, matching SciPy:
// n < 0 is ErrInvalidLength, n = 0 is empty and n = 1 is [1]
// short is true when n is one of these, and w and err are then the window
func shortWindow(n int) (w []float64, short bool, err error) {
	switch {
	case n < 0:
		return nil, true, fmt.Errorf("%w: n must be >=0, got %d", ErrInvalidLength, n)
	case n == 0:
		return []float64{}, true, nil
	case n == 1:
		return []float64{1}, true, nil
	}
	return nil, false, nil
}
