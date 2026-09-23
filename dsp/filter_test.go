package dsp

import (
	"slices"
	"testing"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

func TestFilterRejectsInvalidInput(t *testing.T) {
	cases := []struct {
		name    string
		b, a    []float64
		wantErr string
	}{
		{"empty b", nil, []float64{1}, "b and a must be non-empty"},
		{"empty a", []float64{1}, nil, "b and a must be non-empty"},
		{"a[0] zero", []float64{1}, []float64{0, 1}, "a[0] must be non-zero"},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			_, err := Filter(c.b, c.a, []float64{1, 2, 3})
			assert.ErrorContains(t, err, c.wantErr)
		})
	}
}

// test never mutate inputs
// unnormalised a[0] != 1 case is the one where normalising in place would be tempting
func TestFilterDoesNotModifyInputs(t *testing.T) {
	b, a, x := []float64{2, 1}, []float64{2, -0.5}, []float64{1, 2, 3}
	bBefore, aBefore, xBefore := slices.Clone(b), slices.Clone(a), slices.Clone(x)

	_, err := Filter(b, a, x)
	require.NoError(t, err)

	assert.Equal(t, bBefore, b)
	assert.Equal(t, aBefore, a)
	assert.Equal(t, xBefore, x)
}
