package poly

import (
	"math"
	"testing"

	"github.com/cxd309/tukey/internal/special"
	"github.com/stretchr/testify/assert"
)

// ulp is the distance from |x| to the next float64 up.
func ulp(x float64) float64 {
	x = math.Abs(x)
	return math.Nextafter(x, math.Inf(1)) - x
}

// The big.Float path agrees with float64 where float64 is accurate.
func TestBigDelayMatchesFloat(t *testing.T) {
	g := NewDelayGrid(64)
	c := []float64{1, -0.4, 0.3, 0.1}
	for i := range g.n {
		assert.InDelta(t, g.Delay(c, i), g.bigDelay(c, i), 1e-15, "point %d", i)
	}
}

// bigSinCosPi agrees with special.SinCosPi on every quadrant, within SinCosPi's own few ULP,
// as π·y and math.Sincos each round, where bigSinCosPi rounds once.
func TestBigSinCosPi(t *testing.T) {
	g := NewDelayGrid(7)
	for num := range 4 * g.n {
		bigSin, bigCos := g.bigSinCosPi(num)
		s, _ := bigSin.Float64()
		c, _ := bigCos.Float64()
		wantSin, wantCos := special.SinCosPi(float64(num) / float64(2*g.n))
		assert.InDelta(t, wantSin, s, 5*ulp(wantSin), "sin, num %d", num)
		assert.InDelta(t, wantCos, c, 5*ulp(wantCos), "cos, num %d", num)
	}
}
