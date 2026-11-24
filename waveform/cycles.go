package waveform

import (
	"math"

	"github.com/cxd309/tukey/internal/special"
)

// cosines is cos(2π·(cycles(t) + phi/360)) at each sample time in t,
// for a phase in cycles and an osset phi in degrees.
func cosines(t []float64, phi float64, cycles func(t float64) special.DoubleDouble) (x []float64) {
	offset := special.Exact(phi).Div(special.Exact(360))
	x = make([]float64, len(t))
	for i, ti := range t {
		x[i] = cosCycles(cycles(ti).Add(offset))
	}
	return
}

// cosCycles is cos(2π·c), with c's whole cycles removed first, which is exact,
// so only the fraction left is rounded to a float64.
func cosCycles(c special.DoubleDouble) (cos float64) {
	fraction := (c.Hi - math.Round(c.Hi)) + c.Lo
	_, cos = special.SinCosPi(2 * fraction)
	return
}
