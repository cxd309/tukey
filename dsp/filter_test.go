package dsp

import (
	"slices"
	"testing"

	"github.com/cxd309/godsp/internal/reference"
	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

type filterVector struct {
	reference.Meta
	Params struct {
		B []float64 `json:"b"`
		A []float64 `json:"a"`
		X []float64 `json:"x"`
	} `json:"params"`
	Output struct {
		Y []float64 `json:"y"`
	} `json:"output"`
}

// filterTolerance was measured against scipy.signal.lfilter (see testdata/generate.py):
// FIR and simple-coefficient cases are bit-identical; Butterworth cases agree
// within ~3 ULPs (worst 6.6e-16), differing only in the order of additions.
// Set ~10x above the worst for headroom across platforms.
var filterTolerance = reference.Tolerance{Rel: 1e-14, Abs: 1e-14}

func TestFilterGoldenVectors(t *testing.T) {
	reference.Run(t, "filter", func(t *testing.T, v filterVector) {
		got, err := Filter(v.Params.B, v.Params.A, v.Params.X)
		require.NoError(t, err, v.Description)
		reference.AssertClose(t, "y", got, v.Output.Y, filterTolerance)
	})
}

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
