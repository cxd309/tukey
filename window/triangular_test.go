package window

import (
	"testing"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

func TestTriangularSmall(t *testing.T) {
	w, err := Triang(3)
	require.NoError(t, err)
	assert.Equal(t, []float64{0.5, 1, 0.5}, w)

	w, err = Triang(4)
	require.NoError(t, err)
	assert.Equal(t, []float64{0.25, 0.75, 0.75, 0.25}, w)

	w, err = Bartlett(3)
	require.NoError(t, err)
	assert.Equal(t, []float64{0, 1, 0}, w)

	w, err = Bohman(5)
	require.NoError(t, err)
	assert.Equal(t, 0.0, w[0])
	assert.Equal(t, 0.0, w[4])
}
