package reference

import (
	"encoding/json"
	"math"
	"testing"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

// JSON has no NaN, so the generators write it as null, which Floats reads back as NaN.
func TestFloatsReadsNullAsNaN(t *testing.T) {
	var f Floats
	require.NoError(t, json.Unmarshal([]byte(`[1.5, null, -2, null]`), &f))
	require.Len(t, f, 4)
	assert.Equal(t, 1.5, f[0])
	assert.True(t, math.IsNaN(f[1]))
	assert.Equal(t, -2.0, f[2])
	assert.True(t, math.IsNaN(f[3]))
}

// An empty array is an empty vector, and anything but an array of numbers and nulls is an error.
func TestFloatsDecoding(t *testing.T) {
	var f Floats
	require.NoError(t, json.Unmarshal([]byte(`[]`), &f))
	assert.NotNil(t, f)
	assert.Empty(t, f)
	assert.Error(t, json.Unmarshal([]byte(`["1"]`), &f))
	assert.Error(t, json.Unmarshal([]byte(`{"x": 1}`), &f))
}
