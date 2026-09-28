// Command tolerances reports the error godsp actually achieves
// compares against each SciPy refrence-vector suite
// next to the tolerance its test enforce
//
// the measured numbers are what tolerances and the packaging doc are set from
//
// go run ./cmd/tolerances

package main

import (
	"fmt"
	"os"
	"text/tabwriter"

	"github.com/cxd309/godsp/internal/conformance"
)

// epsilon is the gap between 1.0 and the next float64 (2^-52)
// used to express relative error in ELPs
// it can undercount by up to 2x, sine the relative gap between neighbouring
// flat64s varies between 2^-53 and 2^-52 within each power of two
const epsilon = 0x1p-52

func main() {
	w := tabwriter.NewWriter(os.Stdout, 0, 0, 2, ' ', 0)
	fmt.Fprintln(w, "SUITE\tVECTORS\tFAILING\tWORST\tULPS\tTOLERANCE\tHEADROOM\tWORST CASE")

	failed := false
	for _, suite := range conformance.Suites {
		s, err := summarise(suite)
		if err != nil {
			fmt.Fprintf(os.Stderr, "%s: %v\n", suite.Name, err)
			failed = true
			continue
		}
		failed = failed || s.failing > 0
		fmt.Fprintf(
			w,
			"%s\t%d\t%d\t%.1e\t%.0f\t%.0e\t%s\t%s\n",
			suite.Name,
			s.vectors,
			s.failing,
			s.worst,
			s.worst/epsilon,
			suite.Tolerance.Abs,
			headroom(suite.Tolerance.Abs, s.worst),
			s.worstCase,
		)
	}
	w.Flush()

	if failed {
		os.Exit(1)
	}
}

// summary one suite's measured error across all its refrence vectors
type summary struct {
	vectors   int
	failing   int     // vectors that errored or had an output outside tolerance
	worst     float64 // largest error relative to max|want|, over every output
	worstCase string  // where worst was measured, e.g. "lowpass_order2_wn0.3.json b"
}

// summary run suite and reduce its results to one summary
func summarise(suite conformance.Suite) (s summary, err error) {
	results, err := suite.Run()
	if err != nil {
		return s, err
	}
	s.vectors = len(results)
	for _, r := range results {
		if !r.Passed() {
			s.failing++
		}
		for _, o := range r.Outputs {
			if o.Worst > s.worst {
				s.worst = o.Worst
				s.worstCase = r.File + " " + o.Name
			}
		}
	}
	return
}

// headroom how many times larger the tolerance is than the measured error
func headroom(tolerance, worst float64) (ratio string) {
	if worst == 0 {
		ratio = "exact"
		return
	}
	ratio = fmt.Sprintf("%.0fx", tolerance/worst)
	return
}
