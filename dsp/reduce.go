package dsp

import (
	"fmt"
	"math"
	"slices"

	"github.com/cxd309/tukey/internal/floats"
)

// Sum returns the total of x, or 0 for an empty x.
// It's summed with compensation, so small values beside large ones keep their low bits,
// and it's closer to the exact total than a plain or pairwise sum:
// 1e8, then 20,000 values of 0.001, then -1e8, sum to exactly 20, where MATLAB gives 20.00001.
//
//	total := dsp.Sum(x)
//
// Units: same as x
//
// Equivalent to:
//   - numpy.sum(x), which sums pairwise
//   - MATLAB sum(x), which sums plainly
func Sum(x []float64) (sum float64) {
	var s floats.Sum
	for _, v := range x {
		s.Add(v)
	}
	sum = s.Total()
	return
}

// Mean returns the arithmetic mean of x, its Sum divided by its length,
// summed with compensation, same as Sum.
//
//	m, err := dsp.Mean(x)
//
// Units: same as x
//
// Returns: ErrSignalTooShort for an empty x, which has no mean.
//
// Equivalent to:
//   - numpy.mean(x)
//   - MATLAB mean(x), which gives NaN for an empty x
func Mean(x []float64) (m float64, err error) {
	if len(x) == 0 {
		return 0, fmt.Errorf("%w: Mean of an empty signal", ErrSignalTooShort)
	}
	m = mean(x)
	return
}

// RMS returns the root mean square of x, √(mean(x²)), its squares summed with compensation.
//
//	rms, err := dsp.RMS(x)
//
// Units: same as x
//
// Returns: ErrSignalTooShort for an empty x, which has no RMS.
//
// Equivalent to:
//   - numpy.sqrt(numpy.mean(x**2)), as SciPy has none
//   - MATLAB rms(x), which gives NaN for an empty x
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

// Percentile returns the value below which p percent of x falls, for p from 0 to 100,
// interpolating between the sorted samples, by PercentileLinear unless PercentileHazen is given.
// Any NaN in x gives NaN, as a missing sample leaves the percentile unknown.
// It sorts a copy, so x keeps its order.
//
//	median, err := dsp.Percentile(x, 50)
//
// Units:
//   - p is a percentage
//   - value is in x's unit
//
// Returns:
//   - ErrSignalTooShort for an empty x.
//   - ErrInvalidParameter for a p outside 0 to 100.
//
// Equivalent to:
//   - numpy.percentile(x, p)
//   - MATLAB prctile(x, p), with PercentileHazen, except that prctile ignores NaN
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

// PercentileOption configures Percentile: PercentileLinear or PercentileHazen.
type PercentileOption func(*percentileConfig)

// PercentileLinear interpolates at sorted position p/100·(n-1),
// so 0 and 100 give the minimum and maximum exactly.
// It's the default, and R's and Excel's too.
//
// Equivalent to:
//   - numpy.percentile's method='linear', its default
//   - MATLAB has none
func PercentileLinear() PercentileOption {
	return func(c *percentileConfig) { c.method = percentileLinear }
}

// PercentileHazen interpolates at sorted position p/100·n - 0.5, clamped to the ends.
//
// Equivalent to:
//   - numpy.percentile's method='hazen'
//   - MATLAB prctile
func PercentileHazen() PercentileOption {
	return func(c *percentileConfig) { c.method = percentileHazen }
}

// mean is Sum(x) divided by its length, for use where x is known not to be empty.
// It's NaN for an empty x.
func mean(x []float64) (m float64) {
	m = Sum(x) / float64(len(x))
	return
}

// percentileMethod is where Percentile places p between the sorted samples.
type percentileMethod int

const (
	percentileLinear percentileMethod = iota
	percentileHazen
)

// percentileConfig is Percentile's settings, built from PercentileOptions.
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

// interpolated is sorted's value at a fractional position, from the nearer neighbour, same as numpy's _lerp,
// so it's exact at both neighbours.
// Next to an infinite value it's that infinity, its limit.
func interpolated(sorted []float64, position float64) (value float64) {
	below := math.Floor(position)
	i, fraction := int(below), position-below
	if i+1 >= len(sorted) {
		return sorted[len(sorted)-1]
	}
	a, b := sorted[i], sorted[i+1]
	if v, ok := floats.InfiniteEnd(a, b, fraction); ok {
		return v
	}
	if fraction >= 0.5 {
		return b - (b-a)*(1-fraction)
	}
	return a + (b-a)*fraction
}
