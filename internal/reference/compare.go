package reference

import (
	"fmt"
	"math"
	"math/cmplx"
	"testing"
)

// epsilon is the gap between 1.0 and the next float64 (2^-52)
// used to express relative error in ELPs
// it can undercount by up to 2x, sine the relative gap between neighbouring
// flat64s varies between 2^-53 and 2^-52 within each power of two
const Epsilon = 0x1p-52

// Tolerance says how close got must be to want for a slice comparison
// a value passes if |got-want| <= Abs*max|want| + Rel*|want|
type Tolerance struct {
	Rel float64 // relative to each expected value
	Abs float64 // relative to the largest |want| in the slice, so small but valid slices are still checked
}

// Number is the element types Compare can handle
type Number interface {
	float64 | complex128
}

// Mismatch is one element outside tolerance
type Mismatch struct {
	Index         int
	Diff, Allowed float64
}

// Comparison is the outcome of comparing got against want
type Comparison struct {
	// Worst is the largest |got-want| relative to max|want|:
	// the measured error that tolerances are set from
	Worst      float64
	Mismatches []Mismatch
}

// Compare measures got against want element by element
// NaN in the same posisiton of both counts as a match, as in numpy.allclose(equal_nan=True)
// NaN in only one is always a mismatch
func Compare[T Number](got, want []T, tol Tolerance) (c Comparison, err error) {
	if len(got) != len(want) {
		return c, fmt.Errorf("length mismatch: got %d, want %d", len(got), len(want))
	}
	scale := 0.0
	for _, w := range want {
		if !isNaN(w) {
			scale = math.Max(scale, magnitude(w))
		}
	}
	for i := range want {
		if isNaN(got[i]) && isNaN(want[i]) {
			continue
		}
		diff := magnitude(got[i] - want[i])
		allowed := tol.Abs*scale + tol.Rel*magnitude(want[i])

		relative := diff
		if scale > 0 {
			relative = diff / scale
		}
		c.Worst = max(c.Worst, relative) // the max builtin propagates NaN, so a NaN diff shows here too

		// written as !(within) so a NaN diff is a mismatch
		if !(diff <= allowed) {
			c.Mismatches = append(c.Mismatches, Mismatch{Index: i, Diff: diff, Allowed: allowed})
		}
	}
	return
}

// AssertClose fails t for every element of got outside tol of want
func AssertClose[T Number](t testing.TB, name string, got, want []T, tol Tolerance) {
	t.Helper()
	c, err := Compare(got, want, tol)
	if err != nil {
		t.Fatalf("%s: %v", name, err)
	}
	for _, m := range c.Mismatches {
		t.Errorf("%s[%d]: got %v, want %v (diff %.3e, allowed %.3e)", name, m.Index, got[m.Index], want[m.Index], m.Diff, m.Allowed)
	}
	if testing.Verbose() {
		t.Logf("%s: worst %.3e (%.0f ULPs)", name, c.Worst, c.Worst/Epsilon)
	}
}

// magnitude is |v| for real values and the modulus for complex ones
func magnitude[T Number](v T) (m float64) {
	switch x := any(v).(type) {
	case float64:
		m = math.Abs(x)
	case complex128:
		m = cmplx.Abs(x)
	}
	return
}

// isNaN reports whether v is NaN, or has a NaN part if complex
func isNaN[T Number](v T) (nan bool) {
	switch x := any(v).(type) {
	case float64:
		nan = math.IsNaN(x)
	case complex128:
		nan = cmplx.IsNaN(x)
	}
	return
}
