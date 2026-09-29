// Package conformance checks tukey's public API against the SciPy reference
// vectors in testdata/
//
// Each Suite pairs one dsp function with its vectors and measured tolerance
// the package tests assert on every suite
// cmd/tolerances reports the error each one actually achieves.
package conformance

import (
	"fmt"
	"math"

	"github.com/cxd309/tukey/internal/reference"
)

// Suites is every conformance suite, run by this package's tests and by cmd/tolerances
var Suites = []Suite{
	butterSuite,
	filterSuite,
	filtfiltSuite,
}

// Suite checks one dsp function against one category of reference vectors
type Suite struct {
	Name      string              // the dsp function under test, e.g. "Butter"
	Category  string              // the vectors' directory, testdata/<Category>
	Tolerance reference.Tolerance // measured; see where each suite is defined
	run       func() (results []FileResult, err error)
}

// Output pairs one output of a dsp function with SciPy's for the same input
type Output struct {
	Name      string // e.g. "b"
	Got, Want []float64
	// Scale, if set, is the smallest magnitude errors are measured against;
	// filters set it to their input's, which their rounding error scales with
	// even when the output itself cancels to ~0
	Scale float64
}

// FileResult is how one reference vector file compared
type FileResult struct {
	File    string
	Err     error // the dsp function failed on this vector, or an output's length was wrong
	Outputs []OutputResult
}

// OutputResult is one output's comparison against SciPy
type OutputResult struct {
	Output
	reference.Comparison
}

// Failures describes each element outside tolerance,
// e.g. "b[2]: got 0.1, want 0.2 (diff 1.000e-01, allowed 2.000e-14)"
func (o OutputResult) Failures() (messages []string) {
	for _, m := range o.Mismatches {
		messages = append(messages, fmt.Sprintf("%s[%d]: got %v, want %v (diff %.3e, allowed %.3e)",
			o.Name, m.Index, o.Got[m.Index], o.Want[m.Index], m.Diff, m.Allowed))
	}
	return
}

// Run checks the suite's dsp function against every one of its reference vectors
// err reports vectors that couldn't be loaded; per-vector failures are in results
func (s Suite) Run() (results []FileResult, err error) {
	results, err = s.run()
	return
}

// newSuite builds a Suite from outputs, which runs the dsp function
// on one decoded vector and pairs each of its outputs with SciPy's
func newSuite[T any](name, category string, tol reference.Tolerance, outputs func(v T) ([]Output, error)) (s Suite) {
	s = Suite{Name: name, Category: category, Tolerance: tol}
	s.run = func() (results []FileResult, err error) {
		files, err := reference.Load[T](category)
		if err != nil {
			return nil, err
		}
		for _, f := range files {
			results = append(results, compareFile(f.Name, f.Vector, tol, outputs))
		}
		return
	}
	return
}

// compareFile runs one vector through outputs and compares each result with SciPy's
func compareFile[T any](name string, v T, tol reference.Tolerance, outputs func(v T) ([]Output, error)) (r FileResult) {
	r.File = name
	outs, err := outputs(v)
	if err != nil {
		r.Err = err
		return
	}
	for _, o := range outs {
		c, err := reference.CompareScaled(o.Got, o.Want, tol, o.Scale)
		if err != nil {
			r.Err = fmt.Errorf("%s: %w", o.Name, err)
			return
		}
		r.Outputs = append(r.Outputs, OutputResult{Output: o, Comparison: c})
	}
	return
}

// Passed report whether the vector ran and every output was within tolerance
func (r FileResult) Passed() (passed bool) {
	if r.Err != nil {
		return false
	}
	for _, o := range r.Outputs {
		if len(o.Mismatches) > 0 {
			return false
		}
	}
	return true
}

// maxAbs is the largest magnitude in v
func maxAbs(v []float64) (m float64) {
	for _, x := range v {
		m = math.Max(m, math.Abs(x))
	}
	return
}
