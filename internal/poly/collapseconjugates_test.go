package poly

import (
	"testing"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

// SciPy's _cplxreal docstring example
// pairs collapse to their upper member, sorted by real part then imaginary part, followed by the sorted reals
func TestCollapseConjugates(t *testing.T) {
	roots := []complex128{4, 3, 1, 2 - 2i, 2 + 2i, 2 - 1i, 2 + 1i, 2 - 1i, 2 + 1i, 1 + 1i, 1 - 1i}
	got, err := CollapseConjugates(roots)
	require.NoError(t, err)
	assert.Equal(t, []complex128{1 + 1i, 2 + 1i, 2 + 1i, 2 + 2i, 1, 3, 4}, got)
}
