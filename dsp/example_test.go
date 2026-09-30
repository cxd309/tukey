package dsp_test

import (
	"fmt"
	"log"
	"slices"

	"github.com/cxd309/tukey/dsp"
)

func ExampleButter() {
	// a 2nd-order lowpass at 150Hz, for a signal sampled at 1000 Hz
	f, err := dsp.Butter(2, dsp.Lowpass(150.0/(1000.0/2)))
	if err != nil {
		log.Fatal(err)
	}
	b, a, err := f.BA()
	fmt.Printf("b = %.4f\na = %.4f\n", b, a)
	// Output:
	// b = [0.1311 0.2622 0.1311]
	// a = [1.0000 -0.7478 0.2722]
}

func ExampleFiltFilt() {
	// a narrow pulse at sample 50
	x := make([]float64, 100)
	x[50] = 1

	f, err := dsp.Butter(2, dsp.Lowpass(0.1))
	if err != nil {
		log.Fatal(err)
	}
	b, a, err := f.BA()
	if err != nil {
		log.Fatal(err)
	}
	causal, err := dsp.Filter(b, a, x)
	if err != nil {
		log.Fatal(err)
	}
	zeroPhase, err := dsp.FiltFilt(b, a, x)
	if err != nil {
		log.Fatal(err)
	}

	// Filter delays the pulse; FiltFilt smooths it but leaves it where it was
	fmt.Println("Filter peak at sample", peak(causal))
	fmt.Println("FiltFilt peak at sample", peak(zeroPhase))
	// Output:
	// Filter peak at sample 54
	// FiltFilt peak at sample 50
}

// peak is the index of v's largest value
func peak(v []float64) (index int) {
	index = slices.Index(v, slices.Max(v))
	return
}
