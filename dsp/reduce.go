package dsp

import (
	"fmt"
	"math"
	"slices"

	"github.com/cxd309/tukey/internal/floats"
)

// Sum returns the total of x, 0 for an empty x.
// It's summed with compensation, to maintain low bit accuracy.
//
// Equivalent to numpy.sum(x), which sums pairwise, and MATLAB sum(x), which sums plainly.
// On such signals Sum is closer to the exact total than either:
// one value of 1e8 among 20,000 of 0.001, then cancelled, sum to 20, where MATLAB gives 20.00001
func Sum(x []float64) (sum float64) {
	var s floats.Sum
	for _, v := range x {
		s.Add(v)
	}
	sum = s.Total()
	return
}

// Mean returns the arithmetic mean of x: its Sum divided by its length.
//
// Returns ErrSignalTooShort for an empty x, which has no mean (MATLAB's mean gives NaN)
//
// Equivalent to numpy.mean(x) and MATLAB mean(x), summed with compensation like Sum
func Mean(x []float64) (m float64, err error) {
	if len(x) == 0 {
		return 0, fmt.Errorf("%w: Mean of an empty signal", ErrSignalTooShort)
	}
	m = mean(x)
	return
}

// RMS returns the root mean square of x, sqrt(mean(x^2))
//
// returns ErrSignalTooShort if x is empty (MATLAB's rms gives NaN)
//
// equivalent to MATLAB's rms(x).
// equivalent to numpy.sqrt(numpy.mean(x**2)), no SciPy equivalent
func RMS(x []float64) (rms float64, err error) {
	if len(x) == 0 {
		return 0, fmt.Errorf("%w: RMS of an empty signal", ErrSignalTooShort)
	}
	var sum floats.Sum
	for _, v := range x {
		sum.Add(v * v)
	}
	rms = math.Sqrt(sum.Total() / float64(len(x)))
	return
}

// Percentile returns the value below which p percent of x falls, p from 0 to 100,
// interpolating between the sorted samples, linearly by default, or with PercentileHazen.
//
// Returns:
//   - ErrSignalTooShort for an empty x
//   - ErrInvalidParameter for a p outside 0-100
//
// Equivalent to:
//   - numpy.percentile(x, p)
//   - matlab prctile(x, p) (with PercentileHazen) except that prctile ignores NaN
func Percentile(x []float64, p float64, opts ...PercentileOption) (value float64, err error) {
	var cfg percentileConfig
	for _, opt := range opts {
		opt(&cfg)
	}
	switch {
	case len(x) == 0:
		return 0, fmt.Errorf("%w: Percentile of an empty signal", ErrSignalTooShort)
	case !(p >= 0 && p <= 100):
		return 0, fmt.Errorf("%w: p %g, must be between 0 and 100", ErrInvalidParameter, p)
	}
	if slices.ContainsFunc(x, math.IsNaN) {
		return math.NaN(), nil
	}
	sorted := slices.Clone(x)
	slices.Sort(sorted)
	value = interpolated(sorted, cfg.method.position(p/100, len(sorted)))
	return
}

// PercentileOption configures Percentile
type PercentileOption func(*percentileConfig)

// PercentileLinear interpolates at sorted position p/100·(n-1),
// so 0 and 100 give the minimum and maximum exactly.
// The default, same as numpy, R and Excel.
func PercentileLinear() PercentileOption {
	return func(c *percentileConfig) { c.method = percentileLinear }
}

// PercentileHazen interpolates at sorted position p/100·n - 0.5, clamped to the ends,
// same as matlab prctile and numpy "hazen".
func PercentileHazen() PercentileOption {
	return func(c *percentileConfig) { c.method = percentileHazen }
}

// mean is Sum(x) divided by its length, for internal use where x is known not to be empty.
// NaN for an empty x
func mean(x []float64) (m float64) {
	m = Sum(x) / float64(len(x))
	return
}

// percentileConfdig is Percentile's settings, built from PercentileOptions
// percentileMethod is where Percentile places p between the sorted samples
type percentileMethod int

const (
	percentileLinear percentileMethod = iota
	percentileHazen
)

type percentileConfig struct {
	method percentileMethod
}

// position is where fraction q of the way through n sorted samples falls,
// clamped to the samples, in the same arithmetic order as numpy.
func (m percentileMethod) position(q float64, n int) (position float64) {
	if m == percentileHazen {
		position = float64(n)*q + 0.5 - 1
	} else {
		position = float64(n-1) * q
	}
	position = min(max(position, 0), float64(n-1))
	return
}

// interpolated is sorted's value at a fractional position, from the nearer neighbour,
// sames as numpy lerp does, is it is exact at both neighbours
func interpolated(sorted []float64, position float64) (value float64) {
	below := math.Floor(position)
	i, fraction := int(below), position-below
	if i+1 >= len(sorted) {
		return sorted[len(sorted)-1]
	}
	a, b := sorted[i], sorted[i+1]
	if fraction >= 0.5 {
		return b - (b-a)*(1-fraction)
	}
	return a + (b-a)*fraction
}
