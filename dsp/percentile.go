package dsp

import (
	"fmt"
	"math"
	"slices"
)

// PercentileMethod is how Percentile places p between the sorted samples.
type PercentileMethod int

const (
	// PercentileLinear interpolates at sorted position p/100·(n-1),
	// so 0 and 100 give the minimum ad maximum exactly.
	// The default option, same as numpy, R and Excel
	PercentileLinear PercentileMethod = iota
	// PercentileHazen interpolates at sorted position p/100·n - 0.5,
	// clamped to the ends, same as matlab prctile and numpy "hazen"
	PercentileHazen
)

// PercentileOption configures Percentile
type PercentileOption func(*percentileConfig)

// percentileConfdig is Percentile's settings, built from PercentileOptions
type percentileConfig struct {
	method PercentileMethod
}

// PercentileBy sets how Percentile interpolates.
// PercentileLinear by default.
func PercentileBy(method PercentileMethod) PercentileOption {
	return func(c *percentileConfig) { c.method = method }
}

// Percentile returns the value below which p percent of x falls, p from 0 to 100,
// intpolating between the sorted samples according to PercentileBy.
//
// Returns:
//   - ErrInvalidLength for an empty x
//   - ErrInvalidParameter for a p outside 0-100
//   - ErrInvalidMode for an unknown PercentileMethod
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
		return 0, fmt.Errorf("%w: Percentile of an empty signal", ErrInvalidLength)
	case !(p >= 0 && p <= 100):
		return 0, fmt.Errorf("%w: p %g, must be between 0 and 100", ErrInvalidParameter, p)
	case cfg.method != PercentileLinear && cfg.method != PercentileHazen:
		return 0, fmt.Errorf("%w: PercentileMethod %d", ErrInvalidMode, cfg.method)
	}
	if slices.ContainsFunc(x, math.IsNaN) {
		return math.NaN(), nil
	}
	sorted := slices.Clone(x)
	slices.Sort(sorted)
	value = interpolated(sorted, cfg.method.position(p/100, len(sorted)))
	return
}

// position is where fraction q of the way through n sorted samples falls,
// clamped to the samples, in the same arithmetic order as numpy.
func (m PercentileMethod) position(q float64, n int) (position float64) {
	if m == PercentileHazen {
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
