package fft

import (
	"strconv"
	"testing"

	"github.com/cxd309/tukey/internal/reference"
)

// each dedicated butterfly computes the same p-point DFT as the direct sum
func TestDedicatedButterfliesMatchDirectSum(t *testing.T) {
	dedicated := map[int]func(t, dst []complex128, step int){
		2: butterfly2, 3: butterfly3, 4: butterfly4, 5: butterfly5,
	}
	for p, fast := range dedicated {
		t.Run(strconv.Itoa(p), func(t *testing.T) {
			terms := randomSignal(p)
			want := make([]complex128, p)
			directButterfly(terms, twiddleTable(p), 1, want, 1)

			// step 3, to check the strided writes land in the right places
			got := make([]complex128, 3*p)
			fast(terms, got, 3)
			strided := make([]complex128, p)
			for s := range p {
				strided[s] = got[3*s]
			}
			reference.AssertClose(t, "X", strided, want, reference.Tolerance{Rel: 1e-15, Scaled: 1e-15})
		})
	}
}
