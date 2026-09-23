package dsp

import (
	"encoding/json"
	"os"
	"path/filepath"
	"testing"

	"github.com/stretchr/testify/require"
)

// vectorMeta is the header every golden vector file carries
type vectorMeta struct {
	Description  string `json:"description"`
	SciPyVersion string `json:"scipy_version"`
	NumPyVersion string `json:"numpy_version"`
}

// runGolden decodes every testdata/<category>/*.json into a T
// and runs check on it as a subtest named after the file
func runGolden[T any](t *testing.T, category string, check func(t *testing.T, v T)) {
	t.Helper()
	files, err := filepath.Glob(filepath.Join("testdata", category, "*.json"))
	require.NoError(t, err)
	require.NotEmpty(t, files, "no %s golden vectors; run `just fixtures`", category)

	for _, f := range files {
		t.Run(filepath.Base(f), func(t *testing.T) {
			data, err := os.ReadFile(f)
			require.NoError(t, err)
			var v T
			require.NoError(t, json.Unmarshal(data, &v), "decoding %s", f)
			check(t, v)
		})
	}
}
