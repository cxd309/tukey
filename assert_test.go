package dsp

import (
	"math"
	"math/cmplx"
	"testing"
)

// tolerance describes the difference between got and want
// for a slice comparison a value pass if
// |got-want| <= abs * max|want| + rel*|want|
type tolerance struct {
	rel float64 // relative to each expected value
	abs float64 // relative to the largest |want| in the slice, tiny-but-valid still checked
}

// number is the element types assertAllClose can compare
type number interface {
	float64 | complex128
}

// magnitude is |v| for real values and the modulus for complex ones
func magnitude[T number](v T) (m float64) {
	switch x := any(v).(type) {
	case float64:
		m = math.Abs(x)
	case complex128:
		m = cmplx.Abs(x)
	}
	return
}

// assertAllClose is a numpy.allclose-style comparison, with abs scaled by
// the largest expected magnitude (see tolerance)
// for complex values the distance is the modulus |got-want|
func assertAllClose[T number](t testing.TB, name string, got, want []T, tol tolerance) {
	t.Helper()
	if len(got) != len(want) {
		t.Fatalf("%s: length mismatch: got %d, want %d", name, len(got), len(want))
	}
	scale := 0.0
	for _, w := range want {
		scale = math.Max(scale, magnitude(w))
	}
	for i := range want {
		diff := magnitude(got[i] - want[i])
		allowed := tol.abs*scale + tol.rel*magnitude(want[i])
		if diff > allowed {
			t.Errorf("%s[%d]: got %v, want %v (diff %.3e, allowed %.3e)", name, i, got[i], want[i], diff, allowed)
		}
	}
}
