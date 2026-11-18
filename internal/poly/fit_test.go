package poly

import (
	"math"
	"testing"

	"github.com/cxd309/tukey/internal/floats"
	"github.com/stretchr/testify/assert"
)

// The polynomials stay orthogonal at any degree, up to n - 1:
// q_i·q_j is ‖q_j‖² for i = j and 0 otherwise, to a few ULP.
func TestFitIsOrthogonal(t *testing.T) {
	for _, n := range []int{5, 21, 101} {
		f := NewFit(n, n-1)
		worst := 0.0
		for i := range f.q {
			for j := range f.q {
				want := 0.0
				if i == j {
					want = f.norm[j]
				}
				// relative to ‖q_i‖·‖q_j‖
				worst = math.Max(worst, math.Abs(floats.Dot(f.q[i], f.q[j])-want)/math.Sqrt(f.norm[i]*f.norm[j]))
			}
		}
		assert.Less(t, worst, 1e-14, "n %d", n)
	}
}
