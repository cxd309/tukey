package dsp

import (
	"slices"
	"testing"

	"github.com/stretchr/testify/assert"
)

func TestPaddingExtend(t *testing.T) {
	x := []float64{1, 2, 4, 7}
	cases := []struct {
		name    string
		padding Padding
		n       int
		want    []float64
	}{
		// left: 2*1-4, 2*1-2   right: 2*7-4, 2*7-2
		{"odd", PaddingOdd, 2, []float64{-2, 0, 1, 2, 4, 7, 10, 12}},
		{"even", PaddingEven, 2, []float64{4, 2, 1, 2, 4, 7, 4, 2}},
		{"constant", PaddingConstant, 2, []float64{1, 1, 1, 2, 4, 7, 7, 7}},
		{"longest odd extension, n = len(x)-1", PaddingOdd, 3, []float64{-5, -2, 0, 1, 2, 4, 7, 10, 12, 13}},
		{"none", PaddingNone, 0, []float64{1, 2, 4, 7}},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			xBefore := slices.Clone(x)
			assert.Equal(t, c.want, c.padding.extend(x, c.n))
			assert.Equal(t, xBefore, x, "extend must not modify x")
		})
	}
}

func TestPaddingExtendReturnsNewSlice(t *testing.T) {
	x := []float64{1, 2, 3}
	extended := PaddingNone.extend(x, 0)
	extended[0] = 99
	assert.Equal(t, 1.0, x[0], "extend must not alias x")
}

func TestPaddingExtendPanicsOnBrokenPreconditions(t *testing.T) {
	cases := []struct {
		name    string
		padding Padding
		x       []float64
		n       int
	}{
		{"n = len(x)", PaddingOdd, []float64{1, 2, 3}, 3},
		{"negative n", PaddingOdd, []float64{1, 2, 3}, -1},
		{"empty x", PaddingOdd, nil, 0},
		{"PaddingNone with n > 0", PaddingNone, []float64{1, 2, 3}, 1},
		{"unknown padding", Padding(7), []float64{1, 2, 3}, 1},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			assert.Panics(t, func() { c.padding.extend(c.x, c.n) })
		})
	}
}
