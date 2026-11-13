package dsp

import (
	"math"
	"slices"
	"testing"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

// impulseResponse is SavGolFilter's output for a unit impulse in the middle of
// 13 samples, far enough from the ends that only the centred fit sees it.
// Output i is weight 8 − i of the centred window, for n = 5.
func impulseResponse(t *testing.T, polyOrder int, opts ...SavGolOption) (y []float64) {
	t.Helper()
	x := make([]float64, 13)
	x[6] = 1
	y, err := SavGolFilter(x, 5, polyOrder, opts...)
	require.NoError(t, err)
	return y[4:9]
}

// The textbook coefficients for n = 5, polyOrder 2 (Savitzky and Golay, 1964):
// smoothing (−3, 12, 17, 12, −3)/35, first derivative (−2, −1, 0, 1, 2)/10, and
// second derivative (2, −1, −2, −1, 2)/7, read reversed from the impulse response.
func TestSavGolFilterTextbookCoefficients(t *testing.T) {
	cases := []struct {
		name  string
		deriv int
		want  []float64 // the coefficients, in window order
	}{
		{"smoothing", 0, []float64{-3.0 / 35, 12.0 / 35, 17.0 / 35, 12.0 / 35, -3.0 / 35}},
		{"first derivative", 1, []float64{-0.2, -0.1, 0, 0.1, 0.2}},
		{"second derivative", 2, []float64{2.0 / 7, -1.0 / 7, -2.0 / 7, -1.0 / 7, 2.0 / 7}},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			y := impulseResponse(t, 2, SavGolDeriv(c.deriv, 1))
			assert.InDeltaSlice(t, c.want, reversed(y), 2*ulp(1))
		})
	}
}

// reversed is a copy of v, back to front
func reversed(v []float64) (r []float64) {
	r = slices.Clone(v)
	slices.Reverse(r)
	return
}

// weightSize is the largest Σ|w| of any output's weights, for samples 1 apart,
// which bounds how far the weights' rounding moves an output, relative to x
func weightSize(n, polyOrder, deriv int) (size float64) {
	f := newPolyFit(n, polyOrder)
	for at := range n {
		var sum float64
		for _, w := range f.weights(at, deriv) {
			sum += math.Abs(w)
		}
		size = math.Max(size, sum/math.Pow(float64(f.spacing), float64(deriv)))
	}
	return
}

// A polynomial of degree up to polyOrder is its own fit, at every sample, the
// ends included, and so are its derivatives: here a cubic in t, sampled at
// t = 0, 0.5, 1, ..., whose samples and derivatives are exact.
func TestSavGolFilterReproducesPolynomials(t *testing.T) {
	const dx = 0.5
	p := func(t float64) float64 { return 2*t*t*t - 5*t*t + t - 7 }
	derivatives := []func(t float64) float64{
		p,
		func(t float64) float64 { return 6*t*t - 10*t + 1 },
		func(t float64) float64 { return 12*t - 10 },
		func(t float64) float64 { return 12 },
		func(t float64) float64 { return 0 },
	}
	x := make([]float64, 30)
	largest := 0.0
	for k := range x {
		x[k] = p(float64(k) * dx)
		largest = math.Max(largest, math.Abs(x[k]))
	}
	for _, n := range []int{5, 9, 15} {
		for polyOrder := 3; polyOrder < min(n, 7); polyOrder++ {
			for deriv, want := range derivatives {
				y, err := SavGolFilter(x, n, polyOrder, SavGolDeriv(deriv, dx))
				require.NoError(t, err)
				// a few ULP of the derivative's scale, max |x|/dx^deriv, times
				// the size of its weights, Σ|w|, which for a high derivative
				// at the ends is large: 54 for the 4th of degree 6 at n = 9.
				// Rounding the weights alone costs about half an ULP of that.
				tol := 4 * ulp(1) * largest * weightSize(n, polyOrder, deriv) / math.Pow(dx, float64(deriv))
				for k := range y {
					assert.InDelta(t, want(float64(k)*dx), y[k], tol,
						"n %d, polyOrder %d, deriv %d, y[%d]", n, polyOrder, deriv, k)
				}
			}
		}
	}
}

// A constant is kept exactly, and its derivatives are exactly 0, as each
// output is taken relative to its own sample.
func TestSavGolFilterConstant(t *testing.T) {
	x := make([]float64, 40)
	for k := range x {
		x[k] = 0.1 // not exact in binary: a moving average of it would round
	}
	for _, n := range []int{1, 3, 7, 21} {
		for polyOrder := range min(n, 6) {
			y, err := SavGolFilter(x, n, polyOrder)
			require.NoError(t, err)
			assert.Equal(t, x, y, "n %d, polyOrder %d", n, polyOrder)
			for deriv := 1; deriv <= 3; deriv++ {
				y, err := SavGolFilter(x, n, polyOrder, SavGolDeriv(deriv, 0.01))
				require.NoError(t, err)
				assert.Equal(t, make([]float64, len(x)), y, "n %d, polyOrder %d, deriv %d", n, polyOrder, deriv)
			}
		}
	}
}

// polyOrder 0 is a moving average, and exact: on 1, 2, 3, ..., each centred
// window's mean is its middle sample, and the ends' fit is the mean of the
// first or last n samples.
func TestSavGolFilterMovingAverage(t *testing.T) {
	x := make([]float64, 20)
	for k := range x {
		x[k] = float64(k + 1)
	}
	for _, n := range []int{3, 5, 9} {
		half := n / 2
		want := slices.Clone(x)
		for k := range half {
			want[k] = float64(half + 1)         // the mean of 1 to n
			want[len(x)-1-k] = x[len(x)-1-half] // the mean of the last n
		}
		y, err := SavGolFilter(x, n, 0)
		require.NoError(t, err)
		assert.Equal(t, want, y, "n %d", n)
	}
}

// polyOrder n − 1 fits a polynomial through every sample, so y is x. The
// Gram polynomials' three-term recurrence is off by 1e24 here at n = 101;
// Arnoldi keeps it to a few ULP at any degree.
func TestSavGolFilterThroughEverySample(t *testing.T) {
	x := noise(140, 5)
	for _, n := range []int{3, 5, 21, 101} {
		y, err := SavGolFilter(x, n, n-1)
		require.NoError(t, err)
		assert.InDeltaSlice(t, x, y, 1e-14, "n %d", n)
	}
}

// A small signal on a large offset keeps its digits: the rounding scales with
// the samples' spread about each one, not with the offset, so filtering
// 1e6 + s gives 1e6 + the filtered s, rounded once.
func TestSavGolFilterOffset(t *testing.T) {
	s := noise(100, 6)
	for k := range s {
		s[k] *= 0.01
	}
	offset := make([]float64, len(s))
	for k := range s {
		offset[k] = 1e6 + s[k]
	}
	for _, deriv := range []int{0, 1} {
		ys, err := SavGolFilter(s, 11, 3, SavGolDeriv(deriv, 1))
		require.NoError(t, err)
		y, err := SavGolFilter(offset, 11, 3, SavGolDeriv(deriv, 1))
		require.NoError(t, err)
		for k := range y {
			want := ys[k]
			if deriv == 0 {
				want += 1e6
			}
			// the offset's samples are 1e6 + s rounded, up to half an ULP of
			// 1e6 each, which the fit carries through
			assert.InDelta(t, want, y[k], 2*ulp(1e6), "deriv %d, y[%d]", deriv, k)
		}
	}
}

// Reversing x reverses the fit, and negates its odd derivatives.
func TestSavGolFilterReversal(t *testing.T) {
	x := noise(60, 7)
	for deriv := range 3 {
		y, err := SavGolFilter(x, 9, 4, SavGolDeriv(deriv, 1))
		require.NoError(t, err)
		yr, err := SavGolFilter(reversed(x), 9, 4, SavGolDeriv(deriv, 1))
		require.NoError(t, err)
		sign := 1.0
		if deriv%2 == 1 {
			sign = -1
		}
		for k := range yr {
			yr[k] *= sign
		}
		assert.InDeltaSlice(t, y, reversed(yr), 1e-14, "deriv %d", deriv)
	}
}

// dx scales the deriv-th derivative by 1/dx^deriv.
func TestSavGolFilterSpacing(t *testing.T) {
	x := noise(50, 8)
	for deriv := 1; deriv <= 2; deriv++ {
		unit, err := SavGolFilter(x, 7, 3, SavGolDeriv(deriv, 1))
		require.NoError(t, err)
		y, err := SavGolFilter(x, 7, 3, SavGolDeriv(deriv, 0.25))
		require.NoError(t, err)
		for k := range y {
			assert.InDelta(t, unit[k]*math.Pow(4, float64(deriv)), y[k], 1e-13, "deriv %d, y[%d]", deriv, k)
		}
	}
}

// A derivative of higher order than the fit's degree is 0.
func TestSavGolFilterDerivativeAboveDegree(t *testing.T) {
	y, err := SavGolFilter(noise(20, 9), 5, 2, SavGolDeriv(3, 1))
	require.NoError(t, err)
	assert.Equal(t, make([]float64, 20), y)
}

// With no option, or SavGolDeriv of order 0 at any dx, SavGolFilter smooths.
func TestSavGolFilterDefaultIsSmoothing(t *testing.T) {
	x := noise(40, 10)
	byDefault, err := SavGolFilter(x, 7, 2)
	require.NoError(t, err)
	for _, dx := range []float64{1, 0.001, 50} {
		y, err := SavGolFilter(x, 7, 2, SavGolDeriv(0, dx))
		require.NoError(t, err)
		assert.Equal(t, byDefault, y, "dx %g", dx)
	}
}

// A NaN gives NaN wherever a window holds it: inside, the n outputs centred
// within n/2 of it, and in the first n samples, every output fitted there.
func TestSavGolFilterNaN(t *testing.T) {
	cases := []struct {
		name     string
		at       int
		from, to int // the NaN outputs
	}{
		{"inside", 20, 17, 23},
		{"in the first window", 1, 0, 4},   // the edge fit, 0 to 2, and centred windows to 1 + 3
		{"in the last window", 38, 35, 39}, // centred windows from 38 − 3, and the edge fit
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			x := noise(40, 11)
			x[c.at] = math.NaN()
			y, err := SavGolFilter(x, 7, 2)
			require.NoError(t, err)
			for k, v := range y {
				assert.Equal(t, k >= c.from && k <= c.to, math.IsNaN(v), "y[%d] = %v", k, v)
			}
		})
	}
}

// Infinity sums as it would in plain arithmetic: a lone +Inf makes every
// output that weights it infinite.
func TestSavGolFilterInfinity(t *testing.T) {
	x := []float64{1, 2, 3, 4, math.Inf(1), 6, 7, 8, 9}
	y, err := SavGolFilter(x, 3, 0)
	require.NoError(t, err)
	assert.Equal(t, []float64{2, 2, 3, math.Inf(1), math.Inf(1), math.Inf(1), 7, 8, 8}, y)
}

// n = 1 is the identity, into a new slice.
func TestSavGolFilterOneSampleWindow(t *testing.T) {
	x := noise(10, 12)
	y, err := SavGolFilter(x, 1, 0)
	require.NoError(t, err)
	assert.Equal(t, x, y)
	y[0]++
	assert.NotEqual(t, x[0], y[0], "y shares x's storage")
}

// n = len(x): one window, so every output is the one fit.
func TestSavGolFilterWindowIsWholeSignal(t *testing.T) {
	x := []float64{0, 1, 4, 9, 16} // t²: a quadratic fit reproduces it
	y, err := SavGolFilter(x, 5, 2)
	require.NoError(t, err)
	assert.InDeltaSlice(t, x, y, 4*ulp(16))
}

func TestSavGolFilterLeavesXUnchanged(t *testing.T) {
	x := noise(30, 13)
	original := slices.Clone(x)
	_, err := SavGolFilter(x, 7, 3, SavGolDeriv(1, 0.1))
	require.NoError(t, err)
	assert.Equal(t, original, x)
}

func TestSavGolFilterRejectsInvalidInput(t *testing.T) {
	x := noise(10, 14)
	for _, n := range []int{0, -1, 2, 4} {
		_, err := SavGolFilter(x, n, 0)
		assert.ErrorIs(t, err, ErrInvalidLength, "n %d", n)
	}
	for _, polyOrder := range []int{-1, 5, 6} {
		_, err := SavGolFilter(x, 5, polyOrder)
		assert.ErrorIs(t, err, ErrInvalidParameter, "polyOrder %d", polyOrder)
	}
	_, err := SavGolFilter(x, 5, 2, SavGolDeriv(-1, 1))
	assert.ErrorIs(t, err, ErrInvalidParameter, "deriv −1")
	for _, dx := range []float64{0, -1, math.NaN(), math.Inf(1), math.Inf(-1)} {
		_, err := SavGolFilter(x, 5, 2, SavGolDeriv(1, dx))
		assert.ErrorIs(t, err, ErrInvalidParameter, "dx %g", dx)
	}
	_, err = SavGolFilter(x, 11, 2)
	assert.ErrorIs(t, err, ErrSignalTooShort)
	_, err = SavGolFilter(nil, 1, 0)
	assert.ErrorIs(t, err, ErrSignalTooShort)
}

// The polynomials stay orthogonal at any degree: q_i·q_j is ‖q_j‖² for i = j
// and 0 otherwise, to a few ULP, up to degree n − 1.
func TestPolyFitIsOrthogonal(t *testing.T) {
	for _, n := range []int{5, 21, 101} {
		f := newPolyFit(n, n-1)
		worst := 0.0
		for i := range f.q {
			for j := range f.q {
				want := 0.0
				if i == j {
					want = f.norm[j]
				}
				// relative to ‖q_i‖·‖q_j‖
				worst = math.Max(worst, math.Abs(compensatedDot(f.q[i], f.q[j])-want)/math.Sqrt(f.norm[i]*f.norm[j]))
			}
		}
		assert.Less(t, worst, 1e-14, "n %d", n)
	}
}

// compensatedDot is as if summed exactly and rounded once, where a plain sum
// loses the small term to the large ones, and sums an infinite product as it is.
func TestCompensatedDot(t *testing.T) {
	assert.Equal(t, 1.0, compensatedDot([]float64{1, 1, 1}, []float64{1e16, 1, -1e16}))
	// (1 + 2⁻³⁰)² − 1 = 2⁻²⁹ + 2⁻⁶⁰: the product's rounding error, 2⁻⁶⁰, recovered
	a := 1 + 0x1p-30
	assert.Equal(t, 0x1p-29+0x1p-60, compensatedDot([]float64{a, -1}, []float64{a, 1}))
	assert.Equal(t, math.Inf(1), compensatedDot([]float64{1, 2}, []float64{math.Inf(1), 1}))
	assert.True(t, math.IsNaN(compensatedDot([]float64{1, 1}, []float64{math.Inf(1), math.Inf(-1)})))
}
