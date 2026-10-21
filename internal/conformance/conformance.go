// Package conformance checks tukey's public API against
// reference vectors in testdata/<source>/
//
// Each Suite pairs one function with one source's vectors and measured tolerance;
// a function checked against SciPy and MATLAB has a suite for each,
// sharing the code that runs it
// the package tests assert on every suite
// cmd/tolerances reports the error each one actually achieves.
package conformance

import (
	"fmt"
	"math"

	"github.com/cxd309/tukey/dsp"
	"github.com/cxd309/tukey/internal/reference"
)

// Suite checks one function against one source's category of reference vectors
type Suite struct {
	Name      string              // the function under test, e.g. "Butter"
	Source    reference.Source    // the tool the vectors came from
	Category  string              // the vectors' directory, testdata/<Source>/<Category>
	Tolerance reference.Tolerance // measured; see where each suite is defined
	run       func() (results []FileResult, err error)
}

// Output pairs one output of a tukey function with the reference's for the same input
type Output struct {
	Name      string // e.g. "b"
	Got, Want []float64
	// Scale, if set, is the smallest magnitude errors are measured against;
	// filters set it to their input's, which their rounding error scales with
	// even when the output itself cancels to ~0
	Scale float64
	// GotComplex and WantComplex hold a complex-valued output, in place of Got and Want
	GotComplex, WantComplex []complex128
}

// isComplex reports whether this output holds complex values
func (o Output) isComplex() (ok bool) {
	ok = o.GotComplex != nil || o.WantComplex != nil
	return
}

// compare measures the output against the reference's, real or complex
func (o Output) compare(tol reference.Tolerance) (c reference.Comparison, err error) {
	if o.isComplex() {
		c, err = reference.CompareScaled(o.GotComplex, o.WantComplex, tol, o.Scale)
	} else {
		c, err = reference.CompareScaled(o.Got, o.Want, tol, o.Scale)
	}
	return
}

// values returns element i of got and want, real or complex, for messages
func (o Output) values(i int) (got, want any) {
	if o.isComplex() {
		return o.GotComplex[i], o.WantComplex[i]
	}
	return o.Got[i], o.Want[i]
}

// FileResult is how one reference vector file compared
type FileResult struct {
	File    string
	Err     error // the dsp function failed on this vector, or an output's length was wrong
	Outputs []OutputResult
}

// OutputResult is one output's comparison against the reference
type OutputResult struct {
	Output
	reference.Comparison
}

// Failures describes each element outside tolerance,
// e.g. "b[2]: got 0.1, want 0.2 (diff 1.000e-01, allowed 2.000e-14)"
func (o OutputResult) Failures() (messages []string) {
	for _, m := range o.Mismatches {
		got, want := o.values(m.Index)
		messages = append(messages, fmt.Sprintf("%s[%d]: got %v, want %v (diff %.3e, allowed %.3e)",
			o.Name, m.Index, got, want, m.Diff, m.Allowed))
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
// on one decoded vector and pairs each of its outputs with the vector's
func newSuite[T any](name string, source reference.Source, category string, tol reference.Tolerance, outputs func(v T) ([]Output, error)) (s Suite) {
	s = newSuiteOf(name, source, category, tol, outputs, func(string) bool { return true })
	return
}

// newSuiteOf is newSuite over only the vectors whose file name match accepts.
// It splits a category's hard cases into a suite of their own, so they can have
// a looser tolerance without loosening it for every other vector.
func newSuiteOf[T any](name string, source reference.Source, category string, tol reference.Tolerance, outputs func(v T) ([]Output, error), match func(file string) bool) (s Suite) {
	s = Suite{Name: name, Source: source, Category: category, Tolerance: tol}
	s.run = func() (results []FileResult, err error) {
		files, err := reference.Load[T](source, category)
		if err != nil {
			return nil, err
		}
		for _, f := range files {
			if !match(f.Name) {
				continue
			}
			results = append(results, compareFile(f.Name, f.Vector, tol, outputs))
		}
		return
	}
	return
}

// compareFile runs one vector through outputs and compares each result with the vector's
func compareFile[T any](name string, v T, tol reference.Tolerance, outputs func(v T) ([]Output, error)) (r FileResult) {
	r.File = name
	outs, err := outputs(v)
	if err != nil {
		r.Err = err
		return
	}
	for _, o := range outs {
		c, err := o.compare(tol)
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

// sosOutputs pairs got with want section by section
// each section's b and a as a seperate output
// the overall gain sits only in the first section's b
// this can be many orders of magnitude smaller than the rest,
// so one scale for the whole array would leave it effectively unchecked
func sosOutputs(got dsp.SOS, want [][6]float64) (outputs []Output, err error) {
	if len(got) != len(want) {
		return nil, fmt.Errorf("got %d sections, want %d", len(got), len(want))
	}
	for i := range want {
		outputs = append(outputs,
			Output{Name: fmt.Sprintf("sos[%d].b", i), Got: got[i][:3], Want: want[i][:3]},
			Output{Name: fmt.Sprintf("sos[%d].a", i), Got: got[i][3:], Want: want[i][3:]},
		)
	}
	return
}
