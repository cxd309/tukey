package window_test

import (
	"fmt"
	"log"

	"github.com/cxd309/tukey/window"
)

func ExampleHann() {
	// symmetric, for filter design: zero at both ends
	symmetric, err := window.Hann(5)
	if err != nil {
		log.Fatal(err)
	}
	// periodic, for spectral analysis: the symmetric window one longer, less its last point
	periodic, err := window.Hann(4, window.Periodic())
	if err != nil {
		log.Fatal(err)
	}
	fmt.Println("symmetric:", symmetric)
	fmt.Println("periodic: ", periodic)
	// Output:
	// symmetric: [0 0.5 1 0.5 0]
	// periodic:  [0 0.5 1 0.5]
}
