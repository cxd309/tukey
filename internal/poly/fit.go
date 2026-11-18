package poly

import (
	"math"

	"github.com/cxd309/tukey/internal/floats"
)

// Fit is the polynomials q_0 to q_degree orthogonal over a window of n positions,
// the samples' offsets from its centre divided by its spacing, so from -1 to 1.
// The least-squares fit of a degree is then a projection onto each, with no system to solve:
//
//	Σ q_j·⟨q_j, x⟩/‖q_j‖²
//
// They're built by Arnoldi iteration, same as in Vandermonde with Arnoldi, by Brubeck, Nakatsukasa and Trefethen, 2021.
// Each u·q_j is orthogonalised against every q_i before it, twice,
// which stays orthogonal at any degree.
// The Gram polynomials' three-term recurrence fails beyond a degree of about 2√n.
// q_0 is 1, with ‖q_0‖² = n, so a moving average's weight is 1/n rounded once,
// rather than (1/√n)² rounded three times.
// Every later q_j is normalised, ‖q_j‖² = 1.
type Fit struct {
	n, degree int
	spacing   int         // half the window, or 1 for n = 1
	q         [][]float64 // q[j][k] is q_j at position k
	norm      []float64   // ‖q_j‖², summed over the positions: n, then 1s
	h         [][]float64 // the recurrence, u·q_j = Σ_(i≤j+1) h[i][j]·q_i
}

// NewFit is the orthogonal polynomials up to degree, over n positions.
// degree must be less than n.
func NewFit(n, degree int) (f Fit) {
	f = Fit{n: n, degree: degree, spacing: max(n/2, 1)}
	f.q = make([][]float64, degree+1)
	f.norm = make([]float64, degree+1)
	f.h = make([][]float64, degree+1)
	for j := range f.h {
		f.h[j] = make([]float64, degree)
	}
	f.q[0] = make([]float64, n)
	for k := range f.q[0] {
		f.q[0][k] = 1
	}
	f.norm[0] = float64(n)
	for j := range degree {
		v := make([]float64, n)
		for k := range v {
			v[k] = f.position(k) * f.q[j][k]
		}
		for range 2 {
			for i := 0; i <= j; i++ {
				c := floats.Dot(f.q[i], v) / f.norm[i]
				f.h[i][j] += c
				for k := range v {
					v[k] -= float64(c * f.q[i][k])
				}
			}
		}
		length := math.Sqrt(floats.Dot(v, v))
		for k := range v {
			v[k] /= length
		}
		f.h[j+1][j] = length
		f.q[j+1] = v
		f.norm[j+1] = 1
	}
	return
}

// Spacing is the samples per unit of position:
// half the window, or 1 for a window of 1.
// A derivative in position is spacing^deriv times one in samples.
func (f Fit) Spacing() (spacing int) {
	return f.spacing
}

// Weights maps a window's samples to the fitted polynomial's deriv-th derivative,
// in position u, at sample at:
//
//	Σ_j q_j(u_k)·q_j^(deriv)(u_at)/‖q_j‖²
//
// A derivative above the degree is 0, so every weight is 0.
func (f Fit) Weights(at, deriv int) (w []float64) {
	w = make([]float64, f.n)
	if deriv > f.degree {
		return
	}
	// derivs[m][j] is q_j^(m) at sample at, q_j itself for m = 0,
	// and for m > 0 the recurrence differentiated m times:
	//
	//	h[j+1][j]·q_(j+1)^(m) = u·q_j^(m) + m·q_j^(m-1) - Σ_(i≤j) h[i][j]·q_i^(m)
	u := f.position(at)
	derivs := make([][]float64, deriv+1)
	derivs[0] = make([]float64, f.degree+1)
	for j := range derivs[0] {
		derivs[0][j] = f.q[j][at]
	}
	for m := 1; m <= deriv; m++ {
		derivs[m] = make([]float64, f.degree+1)
		for j := range f.degree {
			var sum floats.Sum
			sum.Add(float64(u * derivs[m][j]))
			sum.Add(float64(float64(m) * derivs[m-1][j]))
			for i := 0; i <= j; i++ {
				sum.Add(-float64(f.h[i][j] * derivs[m][i]))
			}
			derivs[m][j+1] = sum.Total() / f.h[j+1][j]
		}
	}
	coeffs := make([]float64, f.degree+1)
	for j := deriv; j <= f.degree; j++ {
		coeffs[j] = derivs[deriv][j] / f.norm[j]
	}
	for k := range w {
		var sum floats.Sum
		for j := deriv; j <= f.degree; j++ {
			product := float64(f.q[j][k] * coeffs[j])
			sum.Add(product)
			sum.Add(math.FMA(f.q[j][k], coeffs[j], -product))
		}
		w[k] = sum.Total()
	}
	return
}

// position is sample k's position in the window, from -1 to 1.
func (f Fit) position(k int) float64 {
	return float64(k-f.n/2) / float64(f.spacing)
}
