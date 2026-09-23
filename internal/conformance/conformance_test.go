package conformance

import (
	"testing"

	"github.com/stretchr/testify/require"
)

// TestConformance runs every suite, with one subtest per suite and per vector file,
// e.g. TestConformance/Filter/butter2_lp0.3_step.json
func TestConformance(t *testing.T) {
	for _, suite := range Suites {
		t.Run(suite.Name, func(t *testing.T) { checkSuite(t, suite) })
	}
}

// checkSuite runs each of suite's reference vectors as its own subtest
func checkSuite(t *testing.T, suite Suite) {
	results, err := suite.Run()
	require.NoError(t, err)
	for _, r := range results {
		t.Run(r.File, func(t *testing.T) { checkFile(t, r) })
	}
}

// checkFile fails if the dsp function errored on this vector,
// or if any of its outputs fell outside tolerance
func checkFile(t *testing.T, r FileResult) {
	require.NoError(t, r.Err)
	for _, o := range r.Outputs {
		for _, failure := range o.Failures() {
			t.Error(failure)
		}
	}
}
