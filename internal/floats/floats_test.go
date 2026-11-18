package floats

import (
	"math"
	"math/big"
	"testing"

	"github.com/stretchr/testify/assert"
)

func TestSum(t *testing.T) {
	cases := []struct {
		name   string
		values []float64
		want   float64
	}{
		{"empty", nil, 0},
		{"exact", []float64{1, 2, 3}, 6},
		// a plain sum gives 0, as each 1 is lost against 1e100
		{"small lost against large", []float64{1, 1e100, 1, -1e100}, 2},
		// 0.1 isn't exact in binary, so a plain sum of ten drifts to 0.9999999999999999
		{"many small", []float64{0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1}, 1},
		// an infinite value or an overflowing total is infinite, as for a plain sum
		{"infinity", []float64{1, math.Inf(1), 2}, math.Inf(1)},
		{"negative infinity first", []float64{math.Inf(-1), 1}, math.Inf(-1)},
		{"overflow", []float64{math.MaxFloat64, math.MaxFloat64}, math.Inf(1)},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			var s Sum
			for _, v := range c.values {
				s.Add(v)
			}
			assert.Equal(t, c.want, s.Total())
		})
	}
}

// Dot is as if summed exactly and rounded once, where a plain sum loses the small term to the large ones,
// and it sums an infinite product as it is.
func TestDot(t *testing.T) {
	assert.Equal(t, 1.0, Dot([]float64{1, 1, 1}, []float64{1e16, 1, -1e16}))
	// (1 + 2^-30)² - 1 = 2^-29 + 2^-60, so the product's rounding error, 2^-60, is recovered
	a := 1 + 0x1p-30
	assert.Equal(t, 0x1p-29+0x1p-60, Dot([]float64{a, -1}, []float64{a, 1}))
	assert.Equal(t, math.Inf(1), Dot([]float64{1, 2}, []float64{math.Inf(1), 1}))
	assert.True(t, math.IsNaN(Dot([]float64{1, 1}, []float64{math.Inf(1), math.Inf(-1)})))
}

// exactMidpoint is (a + b)/2 for finite a and b, summed exactly in big.Float,
// halved exactly, and rounded once to float64, to nearest even.
func exactMidpoint(a, b float64) float64 {
	sum := new(big.Float).SetPrec(2200).SetFloat64(a)
	sum.Add(sum, new(big.Float).SetFloat64(b))
	half := new(big.Float).SetMantExp(sum, -1)
	m, _ := half.Float64()
	return m
}

// Midpoint is (a + b)/2 rounded once, including where a + b overflows,
// and where the half is subnormal, which a/2 + b/2 rounds twice.
func TestMidpoint(t *testing.T) {
	tiny := math.SmallestNonzeroFloat64
	cases := []struct {
		name string
		a, b float64
		want float64
	}{
		{"exact", 1, 2, 1.5},
		{"overflowing sum", 0x1p1023, 0x1.8p1023, 0x1.4p1023},
		{"overflowing negative sum", -math.MaxFloat64, -math.MaxFloat64, -math.MaxFloat64},
		{"opposite extremes", -math.MaxFloat64, math.MaxFloat64, 0},
		// 1.5 subnormal ULPs, a tie, rounds to even, 2, where a/2 + b/2 gives 0 + 1
		{"subnormal tie", tiny, 2 * tiny, 2 * tiny},
		// 0.5 ULPs, a tie, rounds to even, 0
		{"smallest subnormal and 0", 0, tiny, 0},
		{"infinity", 1, math.Inf(1), math.Inf(1)},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			assert.Equal(t, c.want, Midpoint(c.a, c.b))
		})
	}
	assert.True(t, math.IsNaN(Midpoint(math.Inf(-1), math.Inf(1))), "midpoint of -Inf and +Inf")

	// against the exact midpoint, across magnitudes and signs
	values := []float64{0, tiny, 3 * tiny, 0x1p-1022, 0x1.8p-1022, 0.1, 1, 1 + 0x1p-52, 3, 1e300, math.MaxFloat64}
	for _, a := range values {
		for _, b := range values {
			for _, signs := range [][2]float64{{1, 1}, {1, -1}, {-1, 1}, {-1, -1}} {
				x, y := signs[0]*a, signs[1]*b
				assert.Equal(t, exactMidpoint(x, y), Midpoint(x, y), "Midpoint(%v, %v)", x, y)
			}
		}
	}
}

// Next to an infinite end, interpolation takes the limit, where plain arithmetic gives NaN.
func TestInfiniteEnd(t *testing.T) {
	inf := math.Inf(1)
	cases := []struct {
		name    string
		a, b, t float64
		want    float64
	}{
		{"to +Inf", 1, inf, 0.5, inf},
		{"from -Inf", -inf, 1, 0.5, -inf},
		{"at the finite start", 1, inf, 0, 1},
		{"at the finite end", -inf, 1, 1, 1},
		{"at the infinite start", -inf, 1, 0, -inf},
		{"same infinities", inf, inf, 0.3, inf},
		{"before an infinite start", inf, 1, -2, inf},
		{"beyond a finite end", inf, 1, 3, -inf},
		{"before a finite start", 1, inf, -2, -inf},
		{"beyond an infinite end", 1, inf, 3, inf},
	}
	for _, c := range cases {
		v, ok := InfiniteEnd(c.a, c.b, c.t)
		assert.True(t, ok, c.name)
		assert.Equal(t, c.want, v, c.name)
	}
	v, ok := InfiniteEnd(-inf, inf, 0.5)
	assert.True(t, ok)
	assert.True(t, math.IsNaN(v), "opposite infinities")
	v, _ = InfiniteEnd(inf, inf, 2)
	assert.True(t, math.IsNaN(v), "beyond the same infinities")
	_, ok = InfiniteEnd(1, 2, 0.5)
	assert.False(t, ok, "finite ends")
	_, ok = InfiniteEnd(math.NaN(), 2, 0.5)
	assert.False(t, ok, "NaN isn't infinite")
}
