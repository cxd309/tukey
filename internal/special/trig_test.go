package special

import (
	"math"
	"math/big"
	"testing"

	"github.com/stretchr/testify/assert"
)

// SinCosPi reduces exactly, so whole and half turns are exact however far from 0,
// and quarter turns are √2/2 to an ULP, same as math.Sincos(π/4).
func TestSinCosPi(t *testing.T) {
	cases := []struct{ x, sin, cos float64 }{
		{0, 0, 1},
		{0.5, 1, 0},
		{1, 0, -1},
		{-1.5, 1, 0},
		{1000.25, math.Sqrt2 / 2, math.Sqrt2 / 2},
		{-250.5, -1, 0},
		{1e6 + 0.75, math.Sqrt2 / 2, -math.Sqrt2 / 2},
	}
	for _, c := range cases {
		sin, cos := SinCosPi(c.x)
		assert.InDelta(t, c.sin, sin, 2e-16, "sin(π·%g)", c.x)
		assert.InDelta(t, c.cos, cos, 2e-16, "cos(π·%g)", c.x)
	}
}

// SinCosPiProduct adds back a·b's rounding error.
// The reference reduces the exact product mod 2 in big.Float, where it's exact,
// so only the final rounding remains.
// Without the correction, sin(π·a·b) drifts by 1e-14 here.
func TestSinCosPiProductKeepsDigits(t *testing.T) {
	const a, b = 0.3, 250.5
	exact := new(big.Float).SetPrec(200).Mul(big.NewFloat(a), big.NewFloat(b))
	turns, _ := new(big.Float).Quo(exact, big.NewFloat(2)).Int(nil)
	reduced := exact.Sub(exact, new(big.Float).SetPrec(200).Mul(new(big.Float).SetInt(turns), big.NewFloat(2)))
	r, _ := reduced.Float64()
	wantSin, wantCos := math.Sincos(math.Pi * r)

	sin, cos := SinCosPiProduct(a, b)
	assert.InDelta(t, wantSin, sin, 2e-16)
	assert.InDelta(t, wantCos, cos, 2e-16)
	naive := math.Sin(math.Pi * a * b)
	assert.Greater(t, math.Abs(naive-wantSin), 1e-15, "the naive form should drift here, or the case tests nothing")
}
