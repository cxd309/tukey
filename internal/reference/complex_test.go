package reference

import (
	"encoding/json"
	"testing"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

// Values pairs each real part with the imaginary part at the same index.
func TestComplexValues(t *testing.T) {
	var c Complex
	require.NoError(t, json.Unmarshal([]byte(`{"re": [1, 0, -2], "im": [0, 3, 0.5]}`), &c))
	values, err := c.Values()
	require.NoError(t, err)
	assert.Equal(t, []complex128{1, complex(0, 3), complex(-2, 0.5)}, values)

	empty, err := Complex{}.Values()
	require.NoError(t, err)
	assert.NotNil(t, empty)
	assert.Empty(t, empty)
}

// Parts of different lengths can't be paired.
func TestComplexRejectsUnequalParts(t *testing.T) {
	_, err := Complex{Re: []float64{1, 2}, Im: []float64{0}}.Values()
	assert.Error(t, err)
}
