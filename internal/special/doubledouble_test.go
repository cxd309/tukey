package special

import (
	"math"
	"testing"

	"github.com/stretchr/testify/assert"
)

// relativeError is |got - want|/|want|, from the difference of the two DoubleDoubles.
func relativeError(got, want DoubleDouble) (rel float64) {
	rel = math.Abs(got.Sub(want).Float64() / want.Hi)
	return
}

// SumOf and ProductOf are exact, so Hi + Lo recovers the lost digits.
func TestSumOfAndProductOfAreExact(t *testing.T) {
	s := SumOf(1, 0x1p-60)
	assert.Equal(t, DoubleDouble{Hi: 1, Lo: 0x1p-60}, s)

	p := ProductOf(1+0x1p-30, 1+0x1p-30) // 1 + 2^-29 + 2^-60
	assert.Equal(t, DoubleDouble{Hi: 1 + 0x1p-29, Lo: 0x1p-60}, p)
}

// Div's quotient is the mpmath value to its last digit.
func TestDiv(t *testing.T) {
	third := Exact(1).Div(Exact(3))
	assert.Less(t, relativeError(third, DoubleDouble{0.3333333333333333, 1.850371707708594e-17}), 0x1p-104)
	assert.Equal(t, Exact(1000), Exact(20000).Div(Exact(20)))
}

// The want values are mpmath's at 60 digits, split into Hi and Lo.
func TestExpM1(t *testing.T) {
	cases := []struct {
		x    float64
		want DoubleDouble
	}{
		{1e-20, DoubleDouble{1e-20, 5e-41}},
		{-3e-5, DoubleDouble{-2.9999550004499967e-05, 3.2747415088433923e-22}},
		{0.5, DoubleDouble{0.6487212707001282, -4.731568479435833e-17}},
		{-0.75, DoubleDouble{-0.5276334472589853, 1.7984004434373214e-17}},
		{1, DoubleDouble{1.7182818284590453, -7.747991575210629e-17}},
		{13.1, DoubleDouble{488941.41461545986, 2.8855838860080244e-11}},
		{-20, DoubleDouble{-0.9999999979388464, 4.1293110494709923e-17}},
		{300, DoubleDouble{1.9424263952412558e+130, 1.1143145370175699e+114}},
	}
	for _, c := range cases {
		got := ExpM1(Exact(c.x))
		assert.Less(t, relativeError(got, c.want), 0x1p-94, "ExpM1(%g)", c.x)
	}
}

func TestExpM1Limits(t *testing.T) {
	assert.Equal(t, Exact(0), ExpM1(Exact(0)))
	assert.True(t, math.IsInf(ExpM1(Exact(710)).Hi, 1))
	assert.Equal(t, Exact(-1), ExpM1(Exact(-750)))
	assert.True(t, math.IsNaN(ExpM1(Exact(math.NaN())).Hi))
}

// The want values are mpmath's at 60 digits, split into Hi and Lo.
// Near 1 the logarithm is small, and keeps its relative accuracy.
func TestLog(t *testing.T) {
	cases := []struct {
		x    float64
		want DoubleDouble
	}{
		{1e-300, DoubleDouble{-690.7755278982137, -2.3670096176709832e-14}},
		{0.001, DoubleDouble{-6.907755278982137, -2.1613487097372872e-16}},
		{0.9999999999, DoubleDouble{-1.000000082790371e-10, -4.217583732727628e-27}},
		{1.0000000001, DoubleDouble{1.000000082690371e-10, -4.2169170658954805e-27}},
		{2, DoubleDouble{0.6931471805599453, 2.3190468138462996e-17}},
		{80, DoubleDouble{4.382026634673881, 4.076072943897423e-16}},
		{1e300, DoubleDouble{690.7755278982137, 2.3747660028800243e-14}},
	}
	for _, c := range cases {
		got := Log(Exact(c.x))
		assert.Less(t, relativeError(got, c.want), 0x1p-100, "Log(%g)", c.x)
	}
}

func TestLogLimits(t *testing.T) {
	assert.Equal(t, Exact(0), Log(Exact(1)))
	assert.True(t, math.IsInf(Log(Exact(0)).Hi, -1))
	assert.True(t, math.IsNaN(Log(Exact(-1)).Hi))
	assert.True(t, math.IsInf(Log(Exact(math.Inf(1))).Hi, 1))
}

// Log undoes Exp, so the round trip is within both their errors.
// Only where it's well conditioned: near 0, e^x is 1 plus a part with 32 digits in absolute terms,
// and near the bottom of the range, e^x's Lo is subnormal.
func TestLogOfExp(t *testing.T) {
	for _, x := range []float64{-600, -3.5, -0.25, 0.25, 42, 700} {
		d := SumOf(x, x*0x1p-60)
		assert.Less(t, relativeError(Log(Exp(d)), d), 0x1p-90, "Log(Exp(%g))", x)
	}
}
