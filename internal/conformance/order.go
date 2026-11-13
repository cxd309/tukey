package conformance

import (
	"encoding/json"
	"strings"

	"github.com/cxd309/tukey/dsp"
	"github.com/cxd309/tukey/internal/reference"
)

// orderParams is an order selection's specification, as both sources write it
type orderParams struct {
	Wp            json.RawMessage `json:"wp"`
	Ws            json.RawMessage `json:"ws"`
	RippleDB      float64         `json:"ripple_db"`
	AttenuationDB float64         `json:"attenuation_db"`
}

// orderVector is the exact order selection, with wn matching each band
type orderVector struct {
	reference.Meta
	Params orderParams `json:"params"`
	Output struct {
		Order      int       `json:"order"`
		WnPassband []float64 `json:"wn_passband"`
		WnStopband []float64 `json:"wn_stopband"`
	} `json:"output"`
}

// orderMATLABVector is MATLAB's order selection, with wn in its own convention
type orderMATLABVector struct {
	reference.Meta
	Params orderParams `json:"params"`
	Output struct {
		Order int             `json:"order"`
		Wn    json.RawMessage `json:"wn"`
	} `json:"output"`
}

// orderSelection is ButterOrd, Cheby1Ord, Cheby2Ord or EllipOrd
type orderSelection func(pass, stop dsp.Band, rippleDB, attenuationDB float64, opts ...dsp.OrderOption) (int, dsp.Band, error)

// ButterOrd, Cheby1Ord, Cheby2Ord and EllipOrd, measured on both
// architectures. The SciPy vectors are the exact selection, computed with
// mpmath at 50 digits, rather than scipy.signal's, whose bandstop edges come
// from a numerical minimiser good to about 1e-5; its orders and passband-matched
// edges otherwise agree. Each vector checks the order and both MatchBand
// modes: wn_passband with MatchPassband and wn_stopband with MatchStopband.
// Every order matches, and wn is worst at 7.5e-16 (3 ULP), on EllipOrd's
// highpass_wp0.002_ws0.001_rp1.0_rs30.0 wn_stopband.
//
// Against MATLAB, every order is compared, and matches, and wn where MATLAB
// follows one of tukey's conventions, worst at 3.9e-16:
//   - buttord is stopband-matched, but its bandstop edges come from the same
//     minimiser as SciPy's, so they're compared in a suite of their own,
//     worst at 5.7e-5 on bandstop_wp0.3-0.5_ws0.35-0.45_rp3.0_rs30.0.
//   - cheb1ord and ellipord are passband-matched, except that a bandstop's wn
//     is the passband edges unmoved, which for cheb1ord can miss the stopband
//     specification. Not compared.
//   - cheb2ord returns the stopband edges, which is MatchStopband except for a
//     bandpass, where tukey meets only the tighter edge. Not compared.
var (
	butterOrdSciPySuite = newSuite("ButterOrd", reference.SourceSciPy, "buttord", tolerance(1e-14), orderOutputs(dsp.ButterOrd))
	cheby1OrdSciPySuite = newSuite("Cheby1Ord", reference.SourceSciPy, "cheb1ord", tolerance(1e-14), orderOutputs(dsp.Cheby1Ord))
	cheby2OrdSciPySuite = newSuite("Cheby2Ord", reference.SourceSciPy, "cheb2ord", tolerance(1e-14), orderOutputs(dsp.Cheby2Ord))
	ellipOrdSciPySuite  = newSuite("EllipOrd", reference.SourceSciPy, "ellipord", tolerance(1e-14), orderOutputs(dsp.EllipOrd))

	butterOrdMATLABSuite = newSuiteOf("ButterOrd", reference.SourceMATLAB, "buttord", tolerance(1e-14),
		orderMATLABOutputs(dsp.ButterOrd, dsp.MatchStopband, "lowpass", "highpass", "bandpass"), allFiles)
	butterOrdBandstopMATLABSuite = newSuiteOf("ButterOrd (bandstop)", reference.SourceMATLAB, "buttord", tolerance(1e-3),
		orderMATLABOutputs(dsp.ButterOrd, dsp.MatchStopband, "bandstop"), bandFiles("bandstop"))
	cheby1OrdMATLABSuite = newSuiteOf("Cheby1Ord", reference.SourceMATLAB, "cheb1ord", tolerance(1e-14),
		orderMATLABOutputs(dsp.Cheby1Ord, dsp.MatchPassband, "lowpass", "highpass", "bandpass"), allFiles)
	cheby2OrdMATLABSuite = newSuiteOf("Cheby2Ord", reference.SourceMATLAB, "cheb2ord", tolerance(1e-14),
		orderMATLABOutputs(dsp.Cheby2Ord, dsp.MatchStopband, "lowpass", "highpass", "bandstop"), allFiles)
	ellipOrdMATLABSuite = newSuiteOf("EllipOrd", reference.SourceMATLAB, "ellipord", tolerance(1e-14),
		orderMATLABOutputs(dsp.EllipOrd, dsp.MatchPassband, "lowpass", "highpass", "bandpass"), allFiles)
)

// allFiles accepts every vector
func allFiles(string) bool { return true }

// bandFiles accepts the vectors of one band type, named by it
func bandFiles(btype string) func(file string) bool {
	return func(file string) bool { return strings.HasPrefix(file, btype+"_") }
}

// orderBands is the passband and stopband a specification describes, with
// their band type, inferred from the edges as scipy.signal's *ord do: a
// stopband edge above the passband's is a lowpass, and below it a highpass;
// for a pair, a stopband inside the passband is a bandstop, otherwise a bandpass
func orderBands(p orderParams) (pass, stop dsp.Band, btype string, err error) {
	wp, err := parseWn(p.Wp)
	if err != nil {
		return pass, stop, btype, err
	}
	ws, err := parseWn(p.Ws)
	if err != nil {
		return pass, stop, btype, err
	}
	switch {
	case len(wp) == 1 && wp[0] < ws[0]:
		btype = "lowpass"
	case len(wp) == 1:
		btype = "highpass"
	case wp[0] < ws[0]:
		btype = "bandstop"
	default:
		btype = "bandpass"
	}
	if pass, err = bandFromVector(btype, wp); err != nil {
		return pass, stop, btype, err
	}
	stop, err = bandFromVector(btype, ws)
	return pass, stop, btype, err
}

// orderOutputs pairs a selection's order, and its band in each MatchBand
// mode, with the exact vector's
func orderOutputs(selectOrder orderSelection) func(orderVector) ([]Output, error) {
	return func(v orderVector) (outputs []Output, err error) {
		pass, stop, _, err := orderBands(v.Params)
		if err != nil {
			return nil, err
		}
		modes := []struct {
			name  string
			match dsp.MatchBand
			want  []float64
		}{
			{"wn_passband", dsp.MatchPassband, v.Output.WnPassband},
			{"wn_stopband", dsp.MatchStopband, v.Output.WnStopband},
		}
		for _, m := range modes {
			order, wn, err := selectOrder(pass, stop, v.Params.RippleDB, v.Params.AttenuationDB, dsp.Matching(m.match))
			if err != nil {
				return nil, err
			}
			outputs = append(outputs,
				Output{Name: "order (" + m.name + ")", Got: []float64{float64(order)}, Want: []float64{float64(v.Output.Order)}},
				Output{Name: m.name, Got: wn.Edges(), Want: m.want},
			)
		}
		return
	}
}

// orderMATLABOutputs pairs a selection's order with MATLAB's, and its band,
// selected with match, with MATLAB's wn for the band types listed
func orderMATLABOutputs(selectOrder orderSelection, match dsp.MatchBand, wnCompared ...string) func(orderMATLABVector) ([]Output, error) {
	return func(v orderMATLABVector) (outputs []Output, err error) {
		pass, stop, btype, err := orderBands(v.Params)
		if err != nil {
			return nil, err
		}
		order, wn, err := selectOrder(pass, stop, v.Params.RippleDB, v.Params.AttenuationDB, dsp.Matching(match))
		if err != nil {
			return nil, err
		}
		outputs = []Output{{Name: "order", Got: []float64{float64(order)}, Want: []float64{float64(v.Output.Order)}}}
		for _, compared := range wnCompared {
			if compared != btype {
				continue
			}
			want, err := parseWn(v.Output.Wn)
			if err != nil {
				return nil, err
			}
			outputs = append(outputs, Output{Name: "wn", Got: wn.Edges(), Want: want})
		}
		return
	}
}
