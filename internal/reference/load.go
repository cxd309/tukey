package reference

import (
	"bytes"
	"encoding/json"
	"fmt"
	"os"
	"path/filepath"
	"runtime"
	"testing"
)

// Meta is the header every reference vector file carries: embed it in each vector type
type Meta struct {
	Description  string `json:"description"`
	SciPyVersion string `json:"scipy_version"`
	NumPyVersion string `json:"numpy_version"`
}

// File is one decoded reference vector file
type File[T any] struct {
	Name   string // file name, e.g. "butter2_1p0.3_step.json"
	Vector T
}

// Load decodes every testdata/<category>/*.json into a T, in file-name order
// unknown JSON fields are an error, so a mistyped struct tag fails loudly
// instead of silently decoding as an empty field
func Load[T any](category string) (files []File[T], err error) {
	dir := filepath.Join(testdataDir(), category)
	paths, err := filepath.Glob(filepath.Join(dir, "*.json"))
	if err != nil {
		return nil, err
	}
	if len(paths) == 0 {
		return nil, fmt.Errorf("no reference vectors in %s; run `just fixtures`", dir)
	}
	for _, path := range paths {
		data, err := os.ReadFile(path)
		if err != nil {
			return nil, err
		}
		dec := json.NewDecoder(bytes.NewReader(data))
		dec.DisallowUnknownFields()
		var v T
		if err := dec.Decode(&v); err != nil {
			return nil, fmt.Errorf("decoding %s: %w", path, err)
		}
		files = append(files, File[T]{Name: filepath.Base(path), Vector: v})
	}
	return
}

// Run loads category's reference vectors and runs check on each
// as a subtest named after its file
func Run[T any](t *testing.T, category string, check func(t *testing.T, v T)) {
	t.Helper()
	files, err := Load[T](category)
	if err != nil {
		t.Fatal(err)
	}
	for _, f := range files {
		t.Run(f.Name, func(t *testing.T) { check(t, f.Vector) })
	}
}

// testdataDir is the repository's testdata/, located from this source file rather
// than the working directory, so it's found by both go test (which runs in each
// package's own directory) and go run (which runs wherever it's invoked)
func testdataDir() (dir string) {
	_, thisFile, _, _ := runtime.Caller(0)
	dir = filepath.Join(filepath.Dir(thisFile), "..", "..", "testdata")
	return
}
