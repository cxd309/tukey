<a href="../README.md"><img src="../images/logo.svg" alt="tukey" width="120" align="right"></a>

# Development

How tukey is built and tested, and the conventions its code follows.

## Getting started

Requires Go 1.23 or later, and [just](https://github.com/casey/just) to run the recipes.
The tests need nothing else, as the reference vectors are committed.

```sh
just test         # run every test, on each architecture this machine can run
just lint         # check formatting, go vet, golangci-lint, and the comment style
just fmt          # format everything
just tolerances   # print each function's measured error against SciPy and MATLAB
just bench        # run the FFT benchmarks, six runs each, for benchstat
```

Regenerating the reference vectors needs more:

```sh
just vectors-scipy    # needs uv, which installs the pinned SciPy itself
just vectors-matlab   # needs MATLAB R2025a with the Signal Processing Toolbox
just vectors          # both, SciPy first
```

MATLAB must be on the `PATH`, or set as `MATLAB=/path/to/matlab`.

## How tukey is tested

### Reference vectors

Every public function is checked against SciPy, and against MATLAB where MATLAB has it, on the same inputs.

- `reference/scipy` is a uv project with SciPy 1.18.1 pinned.
  Its `main.py` runs each package's generators, which write `testdata/scipy/<category>/*.json`.
- `reference/matlab` reruns each SciPy vector's inputs through MATLAB's equivalent, and writes `testdata/matlab/<category>/*.json` with the same name.
  A case MATLAB has no equivalent for is skipped, with a comment saying why.
- Where SciPy's own answer is measurably off, the SciPy side is the exact answer instead,
  computed with mpmath at 50 digits or in rational arithmetic, and rounded once.
- Where SciPy has no equivalent, as for `FillMissing`, `Envelope`'s RMS method and `FindSignal`,
  the SciPy side is MATLAB's documented definition, which MATLAB's own vectors then confirm.
- Each file records a description, the tool versions, the inputs and the outputs.
  Complex values are written as `{"re": [...], "im": [...]}`.
- Signals are synthetic, and no real-world data is ever committed.

Regenerating clears each folder first, so a renamed or removed case can't linger.

### Two kinds of test

- **Conformance tests**, in `internal/conformance`, check every public function against the reference vectors, through the public API.
  Each function has a suite per source, listed in `conformance.Suites`.
  `cmd/tolerances` runs the same suites and reports their measured error.
- **Unit tests**, in each package, check what vectors can't:
  invalid input and its errors, that inputs aren't modified, invariants such as FFT round trips and Parseval's theorem,
  and internals that can be checked by hand.

### Tolerances

A value passes if `|got - want| ≤ Scaled·scale + Rel·|want|`,
where scale is the largest expected value, or a floor a suite sets, e.g. a filter's input, whose rounding its output inherits.
NaN matches only NaN.

Each tolerance is set from measurement:
run `just tolerances`, take about 10 times the worst error, rounded to a power of ten,
and record what was measured and why in a comment beside the suite.
[docs/accuracy.md](accuracy.md) is updated from the same numbers.

## Adding a public function

1. Write its cases and SciPy generator, then its MATLAB equivalent if there is one.
   Run `just vectors`, check a value or two by hand, compare the two sources, and commit the vectors.
2. Implement it, following the conventions below.
3. Add its conformance suites to `conformance.Suites`.
4. Write its unit tests: invalid input, input left unmodified, and its invariants.
5. Measure it and set its tolerances.
6. Update its package's `doc.go`, [docs/accuracy.md](accuracy.md), [docs/matlab.md](matlab.md) and [docs/scipy.md](scipy.md),
   and add a runnable example if it shows something users need.

## Conventions

### API

- Packages are flat, holding free functions, named after the SciPy function they match where that reads well.
  A name that would repeat its package reads after it instead, e.g. `fft.Forward`, not `fft.FFT`.
- Every function works on whole signals and is stateless:
  slices in, new slices out, and inputs never modified.
  Block processing is by functions that take and return state, `FilterState` and `SOSFilterState`.
- Each choice is its own option, named after what it chooses, e.g. `dsp.PaddingEven()` or `fft.NormOrtho()`.
  A required choice is a value built by constructors, e.g. `dsp.DetrendLinear()`.
  Options are applied in order, so a later one overrides an earlier one.
- Arguments come in a fixed order:
  - an analysis takes the data first, then required parameters, then options, e.g. `Percentile(x, p)`
  - filtering takes the coefficients first, e.g. `Filter(b, a, x)`, same as SciPy
  - a design takes its order first and its band last, e.g. `Ellip(order, rippleDB, attenuationDB, band)`
- Sample spacing is `dx` for calculus and `fs` for spectral analysis, each positive and finite.
  Design frequencies are fractions of Nyquist, `Freqz` and `GroupDelay` return rad/sample, and spectral frequencies are Hz.
- Return values are named, even single ones, as documentation.

### Empty input and results

- An empty result is always an empty slice, never `nil`.
- A function with one output per input sample gives an empty result for an empty input, e.g. `Diff` or `Filter`.
- A function with no meaningful answer for no samples, or that needs samples to work, returns `ErrSignalTooShort`,
  e.g. `Mean`, `Percentile` or `FiltFilt`.
  `Sum` and `Trapz` of nothing are 0, which is a meaningful answer.

### NaN and infinity

- A NaN or infinite sample is never an error, and affects only the outputs it contributes to, as plain arithmetic would.
- A coefficient of exactly 0 contributes nothing, so 0·NaN never spreads a NaN.
- Interpolating next to an infinite value gives that infinity, its limit.
- A NaN or infinite parameter, coefficient, state or window tap returns an error, unless its doc gives it a meaning.

### Errors and panics

- A caller's mistake returns an error wrapping a sentinel from its package's `errors.go`, so it can be checked with `errors.Is`.
  The same mistake has the same sentinel in every package:

  | Sentinel              | Means                                                            |
  | --------------------- | ---------------------------------------------------------------- |
  | `ErrInvalidLength`    | a length wrong in itself, or lengths that should match and don't |
  | `ErrSignalTooShort`   | too few samples for what's asked, including none                 |
  | `ErrInvalidParameter` | a value outside its domain                                       |
  | `ErrInvalidOption`    | an option that's unknown or doesn't apply                        |

  `dsp` adds its own for bands, orders, coefficients, zeros and poles, sections, filter state, steady state and sample points.
- A bug inside tukey panics, with a clear message.
  Anything a caller can pass in returns an error rather than panicking.
- Range checks are written so NaN fails them, e.g. `!(x > lo && x < hi)`.

### Files

A source file reads in the order an engineer would require the information:

1. the main public function or type, and the types it takes and returns
2. its public options and parameters
3. private types, with their constants and methods
4. private functions

A function and its SOS form share a file, as so much of the code is shared.

Every source file has its own test file, holding the tests of what it defines.
The exceptions are `doc.go`, whose examples are in `example_test.go`, `errors.go`, the conformance harness,
`internal/reference/load.go`, which every conformance run exercises, and `cmd/tolerances`.

Code shared between packages goes in `internal`, e.g. `internal/poly` and `internal/special`.
Code used by one package stays in it, unexported.

### Comments

`just lint` checks the character set and semicolons, in comments and in the markdown docs, and unpaired quotes in comments.
The rest is by review.

A public function's comment has four parts:

```go
// SavGolFilter smooths x with a Savitzky-Golay filter:
// each sample is replaced by a polynomial of degree polyOrder,
// fitted by least squares to the n samples centred on it.
//
//	y, err := dsp.SavGolFilter(x, 21, 3)
//
// Units:
//   - n is a number of samples
//   - y is in x's unit
//
// Returns:
//   - ErrInvalidLength for an n that isn't odd and positive.
//   - ErrSignalTooShort for an x shorter than n.
//
// Equivalent to:
//   - scipy.signal.savgol_filter(x, n, polyOrder)
//   - MATLAB sgolayfilt(x, polyOrder, n)
```

1. A description, with an example where one helps.
2. `Units:`, a bullet per argument or result, or `Units: same as x` where everything is in x's unit.
3. `Returns:`, a bullet per error.
4. `Equivalent to:`, SciPy first, then MATLAB, each with any difference.

A section with one item goes on the heading's line, e.g. `Returns: ErrSignalTooShort for an empty x.`
Each comment stands on its own, so it names its errors rather than pointing at another function's.

Writing rules, for comments and the docs alike:

- At most one sentence per line.
  A long sentence may run over lines, broken at a comma.
- No semicolons and no em dashes.
- A colon ends its line.
- Write "X, same as Y", not "X as Y does".
- Only keyboard characters, or those in `.vscode/maths.code-snippets`.
- No unpaired quote on an indented comment line, which godoc shows as code.
- An equation longer than a short expression goes on its own indented line.

## How tukey was built

tukey was developed with an LLM (Claude). Design decisions, review and most of the code was written directly by me, written from proposal snippets from the LLM.

All functions are directly tested against results from SciPy equivalents, and from MATLAB where it has one, with published tolerances that can be verified with `just tolerances`.
