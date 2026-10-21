package special

import (
	"math"
	"testing"

	"github.com/stretchr/testify/assert"
)

// scipy.special.i0 (SciPy 1.18.1) on each side of the [0, 8] boundary and
// across the large-x expansion, below where SciPy overflows
func TestI0MatchesSciPy(t *testing.T) {
	cases := []struct{ x, want float64 }{
		{0.5, 1.0634833707413236},
		{1, 1.2660658777520082},
		{7.999, 427.16443131520265},
		{8, 427.56411572180474},
		{8.001, 427.9641777084088},
		{14, 129418.56270064856},
		{30, 781672297823.9775},
		{100, 1.0737517071310738e+42},
		{709, 1.2315477067016538e+306},
	}
	for _, c := range cases {
		assert.InEpsilon(t, c.want, I0(c.x), 1e-15, "I0(%g)", c.x)
	}
}

// the defining series, the sum of ((x/2)^k / k!)^2, is an independent check:
// every term is positive, so summing it loses nothing to cancellation
func TestI0MatchesSeries(t *testing.T) {
	for x := 0.0; x <= 30; x += 0.25 {
		assert.InEpsilon(t, series(x), I0(x), 1e-14, "I0(%g)", x)
	}
}

// near and beyond exp's overflow (~709.44 on amd64, ~709.78 elsewhere), I0 is
// finite up to ~713.99; checked against the asymptotic expansion
// e^x / sqrt(2 pi x) (1 + 1/(8x) + 9/(128x^2)), whose next term is ~2e-10 here,
// compared in logs so neither side overflows
func TestI0NearExpOverflow(t *testing.T) {
	for x := 700.0; x <= 713.9; x += 0.1 {
		want := x - 0.5*math.Log(2*math.Pi*x) + math.Log1p(1/(8*x)+9/(128*x*x))
		assert.InDelta(t, want, math.Log(I0(x)), 1e-9, "I0(%g)", x)
	}
	assert.True(t, math.IsInf(I0(714), 1), "I0(714) overflows")
}

func TestI0Special(t *testing.T) {
	assert.InEpsilon(t, 1.0, I0(0), 1e-16)
	assert.True(t, math.IsInf(I0(math.Inf(1)), 1))
	assert.True(t, math.IsInf(I0(math.Inf(-1)), 1))
	assert.True(t, math.IsNaN(I0(math.NaN())))
}

// I0 is even: it works on |x|, so exactly
func TestI0IsEven(t *testing.T) {
	for _, x := range []float64{0.5, 8, 14, 100, 712} {
		assert.Equal(t, I0(x), I0(-x), "I0(%g)", x)
	}
}

// series is I0 by its definition, for small enough x that the terms stay finite
func series(x float64) (sum float64) {
	q := x * x / 4
	term := 1.0
	for k := 1; term > sum*1e-17 || k < 2; k++ {
		sum += term
		term *= q / float64(k*k)
	}
	return
}
