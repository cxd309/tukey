package conformance

import (
	"fmt"
	"strings"

	"github.com/cxd309/tukey/internal/reference"
	"github.com/cxd309/tukey/waveform"
)

type chirpVector struct {
	reference.Meta
	Params struct {
		T          []float64 `json:"t"`
		F0         float64   `json:"f0"`
		T1         float64   `json:"t1"`
		F1         float64   `json:"f1"`
		Method     string    `json:"method"`
		Phi        float64   `json:"phi"`
		VertexZero bool      `json:"vertex_zero"`
	} `json:"params"`
	Output struct {
		X []float64 `json:"x"`
	} `json:"output"`
}

// Chirp, measured on both architectures.
// The SciPy vectors are the exact chirp, computed with mpmath at 50 digits and rounded once, rather than scipy.signal's output.
// tukey is within 1 ULP of them on every vector, where scipy.signal is up to 1.8e-10 away on the audio grid,
// as it takes the cosine of the phase in radians, which loses digits as the phase grows.
//
// MATLAB's chirp does the same, and its vectors skip the sweeps it doesn't have:
// hyperbolic, the quadratic with its vertex at t1, and the logarithmic at f0 = f1.
// Over the unit grid it's worst at 8.0e-14, on linear_unit_down.
// Where the phase runs to many cycles, on the beyond and audio grids, it's worst at 1.3e-10, on linear_audio,
// so those vectors have a looser suite of their own.
// Tolerance: about 10 times the worst.
var (
	chirpSciPySuite  = newSuite("Chirp", reference.SourceSciPy, "waveform/chirp", tolerance(2e-15), chirpOutputs)
	chirpMATLABSuite = newSuiteOf("Chirp", reference.SourceMATLAB, "waveform/chirp", tolerance(1e-12), chirpOutputs,
		func(file string) bool { return !isLongChirp(file) })
	chirpLongMATLABSuite = newSuiteOf("Chirp (long)", reference.SourceMATLAB, "waveform/chirp", tolerance(2e-9),
		chirpOutputs, isLongChirp)
)

// isLongChirp reports whether a vector's phase runs to many cycles:
// the beyond grid, out to 3·t1, and the audio grid, ten seconds from 20 Hz to 20 kHz.
func isLongChirp(file string) (ok bool) {
	ok = strings.Contains(file, "_beyond") || strings.Contains(file, "_audio")
	return
}

// chirpMethods is each SciPy method's Option.
var chirpMethods = map[string]func() waveform.Option{
	"linear":      waveform.Linear,
	"quadratic":   waveform.Quadratic,
	"logarithmic": waveform.Logarithmic,
	"hyperbolic":  waveform.Hyperbolic,
}

// chirpOutputs pairs Chirp's output with a vector's.
func chirpOutputs(v chirpVector) (outputs []Output, err error) {
	method, ok := chirpMethods[v.Params.Method]
	if !ok {
		return nil, fmt.Errorf("unknown chirp method %q", v.Params.Method)
	}
	opts := []waveform.Option{method(), waveform.Phase(v.Params.Phi)}
	if !v.Params.VertexZero {
		opts = append(opts, waveform.VertexAtT1())
	}
	x, err := waveform.Chirp(v.Params.T, v.Params.F0, v.Params.T1, v.Params.F1, opts...)
	if err != nil {
		return nil, err
	}
	// Scale 1, as a chirp swings between -1 and 1, and its zero crossings would otherwise have to match exactly
	outputs = []Output{{Name: "x", Got: x, Want: v.Output.X, Scale: 1}}
	return
}

type sweepPolyVector struct {
	reference.Meta
	Params struct {
		T    []float64 `json:"t"`
		Poly []float64 `json:"poly"`
		Phi  float64   `json:"phi"`
	} `json:"params"`
	Output struct {
		X []float64 `json:"x"`
	} `json:"output"`
}

// SweepPoly, measured on both architectures.
// The SciPy vectors are the exact sweep, computed with mpmath at 50 digits and rounded once, rather than scipy.signal's output,
// tukey is within 1 ULP of them on every vector, worst on audio,
// where scipy.signal is up to 1.2e-10 away, for the same reason as Chirp.
// Tolerance: about 10 times the worst.
// MATLAB has none.
var sweepPolySciPySuite = newSuite("SweepPoly", reference.SourceSciPy, "waveform/sweep_poly", tolerance(2e-15), sweepPolyOutputs)

// sweepPolyOutputs pairs SweepPoly's output with a vector's.
func sweepPolyOutputs(v sweepPolyVector) (outputs []Output, err error) {
	x, err := waveform.SweepPoly(v.Params.T, v.Params.Poly, waveform.Phase(v.Params.Phi))
	if err != nil {
		return nil, err
	}
	// Scale 1, same as Chirp
	outputs = []Output{{Name: "x", Got: x, Want: v.Output.X, Scale: 1}}
	return
}

type squareVector struct {
	reference.Meta
	Params struct {
		T    []float64 `json:"t"`
		Duty float64   `json:"duty"`
	} `json:"params"`
	Output struct {
		X []float64 `json:"x"`
	} `json:"output"`
}

// Square, measured on both architectures.
// The SciPy vectors are scipy.signal's own output, which reduces t by the float64 2π, same as tukey,
// so every sample matches exactly.
//
// MATLAB's mod puts a t within round-off of 2π·k at the start of a cycle,
// so it's +1 where SciPy and tukey are -1 at the 2 samples on the steps grid
// one ULP beyond -4π and -2π, at every duty but 0 and 1.
// Those vectors are left out of its suite, and it matches exactly on the rest.
var (
	squareSciPySuite  = newSuite("Square", reference.SourceSciPy, "waveform/square", reference.Tolerance{}, squareOutputs)
	squareMATLABSuite = newSuiteOf("Square", reference.SourceMATLAB, "waveform/square", reference.Tolerance{}, squareOutputs,
		func(file string) bool { return !strings.HasPrefix(file, "steps_") })
)

// squareOutputs pairs Square's output with a vector's.
func squareOutputs(v squareVector) (outputs []Output, err error) {
	x, err := waveform.Square(v.Params.T, waveform.Duty(v.Params.Duty))
	if err != nil {
		return nil, err
	}
	outputs = []Output{{Name: "x", Got: x, Want: v.Output.X, Scale: 1}}
	return
}

type sawtoothVector struct {
	reference.Meta
	Params struct {
		T    []float64 `json:"t"`
		Duty float64   `json:"duty"`
	} `json:"params"`
	Output struct {
		X []float64 `json:"x"`
	} `json:"output"`
}

// Sawtooth, measured on both architectures.
// The SciPy vectors are exact, from mpmath, reducing t by the float64 2π, same as tukey,
// and tukey is within an ULP of them.
//
// MATLAB's rem puts a t within round-off of 2π·k at the start of a cycle,
// so at duty 0 and 1, where the wave jumps by 2 at the start of a cycle,
// it's on the other side of the jump at the 2 samples on the steps grid one ULP beyond -4π and -2π.
// Those vectors are left out of its suite.
// On the rest MATLAB is 9.0e-15 from tukey, from scaling t to a fraction of a cycle before its ramps.
var (
	sawtoothSciPySuite  = newSuite("Sawtooth", reference.SourceSciPy, "waveform/sawtooth", tolerance(2e-15), sawtoothOutputs)
	sawtoothMATLABSuite = newSuiteOf("Sawtooth", reference.SourceMATLAB, "waveform/sawtooth", tolerance(1e-13), sawtoothOutputs,
		func(file string) bool { return file != "steps_duty0.json" && file != "steps_duty1.json" })
)

// sawtoothOutputs pairs Sawtooth's output with a vector's.
func sawtoothOutputs(v sawtoothVector) (outputs []Output, err error) {
	x, err := waveform.Sawtooth(v.Params.T, waveform.Duty(v.Params.Duty))
	if err != nil {
		return nil, err
	}
	outputs = []Output{{Name: "x", Got: x, Want: v.Output.X, Scale: 1}}
	return
}
