package dsp

import (
	"math"
	"math/cmplx"
	"testing"

	"github.com/stretchr/testify/assert"
)

// designPrototypes are a prototype without finite zeros, and one with them, so every zero rule is exercised.
var designPrototypes = []struct {
	name  string
	proto analogZPK
}{
	{"butterworth", butterworthPrototype(3)},
	{"chebyshev2", chebyshev2Prototype(4, 40)},
}

// analogResponse is f's response at s, Gain·Π(s - zero)/Π(s - pole).
func analogResponse(f analogZPK, s complex128) (h complex128) {
	h = complex(f.Gain, 0)
	for _, z := range f.Zeros {
		h *= s - z
	}
	for _, p := range f.Poles {
		h /= s - p
	}
	return
}

// digitalResponse is f's response at z, Gain·Π(z - zero)/Π(z - pole).
func digitalResponse(f ZPK, z complex128) (h complex128) {
	return analogResponse(analogZPK(f), z)
}

// assertSameResponse checks got against want relative to want's magnitude.
func assertSameResponse(t *testing.T, want, got complex128, msg string, args ...any) {
	t.Helper()
	assert.InDelta(t, 0, cmplx.Abs(got-want)/cmplx.Abs(want), 1e-12, append([]any{msg}, args...)...)
}

// Each frequency transform is a substitution for s,
// so the transformed filter's response at s is the prototype's at the substituted point.
func TestFrequencyTransformsSubstitute(t *testing.T) {
	const wo, bw = 3.0, 0.5
	for _, c := range designPrototypes {
		transforms := []struct {
			name       string
			analog     analogZPK
			substitute func(s complex128) complex128
		}{
			{"lowpass", c.proto.toLowpass(wo), func(s complex128) complex128 { return s / wo }},
			{"highpass", c.proto.toHighpass(wo), func(s complex128) complex128 { return wo / s }},
			{"bandpass", c.proto.toBandpass(wo, bw), func(s complex128) complex128 { return (s*s + wo*wo) / (s * bw) }},
			{"bandstop", c.proto.toBandstop(wo, bw), func(s complex128) complex128 { return s * bw / (s*s + wo*wo) }},
		}
		for _, tr := range transforms {
			for _, w := range []float64{0.3, 1, 2.9, 7} {
				s := complex(0, w)
				want := analogResponse(c.proto, tr.substitute(s))
				assertSameResponse(t, want, analogResponse(tr.analog, s), "%s %s at ω = %g", c.name, tr.name, w)
			}
		}
	}
}

// Each transform keeps the filter's degree, and a bandpass or bandstop doubles it.
func TestFrequencyTransformsKeepDegree(t *testing.T) {
	for _, c := range designPrototypes {
		n := len(c.proto.Poles)
		assert.Len(t, c.proto.toLowpass(2).Poles, n, c.name)
		assert.Len(t, c.proto.toHighpass(2).Zeros, n, c.name)
		assert.Len(t, c.proto.toBandpass(2, 1).Poles, 2*n, c.name)
		assert.Len(t, c.proto.toBandstop(2, 1).Zeros, 2*n, c.name)
	}
}

// The bilinear transform substitutes s = fs2·(z - 1)/(z + 1),
// so the digital response on the unit circle is the analog one at the warped frequency,
// and each zero at infinity lands at z = -1, Nyquist.
func TestBilinearSubstitutes(t *testing.T) {
	for _, c := range designPrototypes {
		digital := c.proto.toLowpass(2).bilinear()
		assert.Len(t, digital.Zeros, len(digital.Poles), c.name)
		for _, theta := range []float64{0.1, 1, 2.5} {
			z := cmplx.Exp(complex(0, theta))
			s := fs2 * (z - 1) / (z + 1)
			assertSameResponse(t, analogResponse(c.proto.toLowpass(2), s), digitalResponse(digital, z), "%s at θ = %g", c.name, theta)
		}
	}
	nyquistZeros := 0
	for _, z := range butterworthPrototype(3).bilinear().Zeros {
		if z == -1 {
			nyquistZeros++
		}
	}
	assert.Equal(t, 3, nyquistZeros, "a Butterworth's three zeros at infinity")
}

// prewarp undoes the bilinear transform's warping,
// so a digital edge at wn maps to the analog frequency the transform sends to wn.
func TestPrewarp(t *testing.T) {
	for _, wn := range []float64{0.01, 0.3, 0.9} {
		assert.InEpsilon(t, 4*math.Tan(math.Pi*wn/2), prewarp(wn), 1e-15, "wn %g", wn)
		z := cmplx.Exp(complex(0, math.Pi*wn))
		assert.InEpsilon(t, prewarp(wn), imag(fs2*(z-1)/(z+1)), 1e-13, "wn %g", wn)
	}
}

// Each split pair r ± √(r² - wo²) is the two roots of x² - 2r·x + wo², so they sum to 2r and multiply to wo².
func TestSplitRoots(t *testing.T) {
	roots := []complex128{complex(-0.5, 0.8), complex(-1, 0), complex(-0.1, -2)}
	const wo = 1.5
	split := splitRoots(roots, wo)
	for i, r := range roots {
		upper, lower := split[i], split[len(roots)+i]
		assertSameResponse(t, 2*r, upper+lower, "sum of root %d's pair", i)
		assertSameResponse(t, complex(wo*wo, 0), upper*lower, "product of root %d's pair", i)
	}
}

// A prototype's poles sit at the angles prototypeAngle gives, symmetric about the real axis.
func TestPrototypeAngleIsSymmetric(t *testing.T) {
	for _, order := range []int{1, 2, 5, 8} {
		for k := range order {
			assert.InDelta(t, -prototypeAngle(order, order-1-k), prototypeAngle(order, k), 1e-15, "order %d, k %d", order, k)
		}
	}
	assert.Equal(t, 0.0, prototypeAngle(3, 1), "an odd order's middle pole is real")
}
