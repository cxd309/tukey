package dsp

import (
	"encoding/json"
	"fmt"
	"os"
	"path/filepath"
	"testing"

	"github.com/stretchr/testify/require"
)

type butterVector struct {
	Description  string `json:"description"`
	SciPyVersion string `json:"scipy_version"`
	NumPyVersion string `json:"numpy_version"`
	Params       struct {
		Order int             `json:"order"`
		Wn    json.RawMessage `json:"wn"`
		BType string          `json:"btype"`
	} `json:"params"`
	Output struct {
		B []float64 `json:"b"`
		A []float64 `json:"a"`
	} `json:"output"`
}

// parseWn handles the golden vectors' wn field, which is a scalar for
// low/high-pass and a [low, high] pair for band-pass/stop.
func parseWn(raw json.RawMessage) (wn []float64, err error) {
	var scalar float64
	if err := json.Unmarshal(raw, &scalar); err == nil {
		return []float64{scalar}, nil
	}
	var pair []float64
	if err := json.Unmarshal(raw, &pair); err == nil {
		return pair, nil
	}
	return nil, fmt.Errorf("wn is neither a float nor a [low, high] pair: %s", raw)
}

func btypeFromString(s string) (bt BandType, err error) {
	switch s {
	case "lowpass":
		return LowPass, nil
	case "highpass":
		return HighPass, nil
	case "bandpass":
		return BandPass, nil
	case "bandstop":
		return BandStop, nil
	default:
		return 0, fmt.Errorf("unknown btype %q", s)
	}
}

var butterTolerance = tolerance{rel: 1e-9, abs: 1e-9}

func TestButterGoldenVectors(t *testing.T) {
	files, err := filepath.Glob("testdata/butter/*.json")
	require.NoError(t, err)
	require.NotEmpty(t, files, "no golden vectors found — run `uv run --project testdata generate.py`")

	for _, f := range files {
		f := f
		t.Run(filepath.Base(f), func(t *testing.T) {
			data, err := os.ReadFile(f)
			require.NoError(t, err)

			var v butterVector
			require.NoError(t, json.Unmarshal(data, &v), "decoding %s", f)

			wn, err := parseWn(v.Params.Wn)
			require.NoError(t, err)
			bt, err := btypeFromString(v.Params.BType)
			require.NoError(t, err)

			var gotB, gotA []float64
			if len(wn) == 1 {
				gotB, gotA, err = Butter(v.Params.Order, wn[0], bt)
			} else {
				gotB, gotA, err = ButterBand(v.Params.Order, wn[0], wn[1], bt)
			}
			require.NoError(t, err, v.Description)
			require.Len(t, gotB, len(v.Output.B), "b length mismatch")
			require.Len(t, gotA, len(v.Output.A), "a length mismatch")

			assertAllClose(t, "b", gotB, v.Output.B, butterTolerance)
			assertAllClose(t, "a", gotA, v.Output.A, butterTolerance)
		})
	}
}
