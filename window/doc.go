// Package window computes window functions, matching scipy.signal.windows.
// It has every SciPy window except dpss, and uses only the standard library.
//
// # Windows
//
// Each window returns n samples:
// symmetric by default, for filter design,
// or periodic, DFT-even, with Periodic(), for spectral analysis.
//
//	Cosine sums:   Hann, Hamming, Blackman, BlackmanHarris, Nuttall, FlatTop,
//	               and the general forms GeneralHamming and GeneralCosine
//	Piecewise:     Boxcar, Triang, Bartlett, Parzen, Bohman, Barthann, Tukey
//	Sine:          Cosine, Lanczos
//	Exponential:   Gaussian, GeneralGaussian, Exponential
//	Bessel:        Kaiser, KaiserBesselDerived, the MDCT window of AAC and Vorbis
//	Low sidelobe:  Taylor, Chebwin, the Dolph-Chebyshev window
//
// # Options
//
// One Option type serves every window.
// A window given an option that doesn't apply to it returns ErrInvalidOption.
//   - Periodic() gives the symmetric window of n+1 points without its last point, for every window but KaiserBesselDerived.
//     It's SciPy's sym=False and MATLAB's 'periodic'.
//   - Center(c) moves Exponential's peak to sample c, with Periodic only.
//   - Unnormalised() leaves Taylor's peak where its formula puts it, rather than at 1.
//     It's SciPy's norm=False, and MATLAB's taylorwin.
//
// # Compatibility
//
// Every window matches its scipy.signal.windows counterpart, and its MATLAB
// R2025a counterpart, within a measured tolerance,
// checked on identical inputs at n = 0, 1, 2, 3, 8, 9, 64 and 255:
//
//	Function             SciPy                  MATLAB          SciPy    MATLAB   Tested to
//	Hann                 hann                   hann            4.4e-16  1.7e-16  1e-14
//	Hamming              hamming                hamming         4.4e-16  2.2e-16  1e-14
//	Blackman             blackman               blackman        5.6e-16  3.3e-16  1e-14
//	BlackmanHarris       blackmanharris         blackmanharris  6.7e-16  7.8e-16  1e-14
//	Nuttall              nuttall                nuttallwin      6.7e-16  8.9e-16  1e-14
//	FlatTop              flattop                flattopwin      7.5e-16  3.3e-16  1e-14
//	GeneralHamming       general_hamming                        4.4e-16           1e-14
//	GeneralCosine        general_cosine                         1.7e-15           1e-13
//	Boxcar               boxcar                 rectwin         exact    exact    exact
//	Triang               triang                 triang          exact    exact    exact
//	Bartlett             bartlett               bartlett        1.1e-16  exact    1e-14
//	Parzen               parzen                 parzenwin       2.2e-16  1.1e-16  1e-14
//	Bohman               bohman                 bohmanwin       5.6e-16  2.2e-16  1e-14
//	Barthann             barthann               barthannwin     3.3e-16  3.3e-16  1e-14
//	Tukey                tukey                  tukeywin        6.7e-16  6.7e-16  1e-14
//	Cosine               cosine                                 4.6e-16           1e-14
//	Lanczos              lanczos                                2.2e-16           1e-14
//	Gaussian             gaussian               gausswin        1.1e-16  1.7e-16  1e-14
//	GeneralGaussian      general_gaussian                       1.1e-16           1e-14
//	Exponential          exponential                            1.1e-16           1e-14
//	Kaiser               kaiser                 kaiser          1.7e-15  2.0e-15  1e-13
//	KaiserBesselDerived  kaiser_bessel_derived                  exact             1e-13
//	Taylor               taylor                 taylorwin       5.6e-16  3.4e-16  1e-14
//	Chebwin              chebwin                chebwin         1.5e-12  3.7e-13  1e-11
//
// Errors are the largest difference, relative to 1 or to the window's peak if that's larger.
// 1e-16 is about 1 ULP of a float64.
// MATLAB's rectwin, triang, bartlett, parzenwin, bohmanwin, barthannwin,
// tukeywin, gausswin, kaiser, taylorwin and chebwin have no periodic form,
// so they're checked symmetric only.
//
// Kaiser's error grows with beta, as I0's sensitivity to its argument does.
// Chebwin's comes from evaluating a Chebyshev polynomial of degree n-1 near ±1,
// where an ULP of difference in its argument is magnified.
// SciPy and MATLAB differ from each other by as much.
//
// # Differences from SciPy and MATLAB
//
// Symmetric windows are computed for their first half and mirrored, so they're exactly symmetric.
// SciPy's are symmetric only to within a few ULP.
//
// MATLAB's hanning is Hann without its zero ends, Hann(n+2) less its first and last sample.
// Hann's doc shows how to build it.
//
// Taylor takes its sidelobe level as positive dB, same as SciPy.
// MATLAB's taylorwin takes it negative and doesn't normalise,
// so taylorwin(n, nbar, -sll) is Taylor(n, nbar, sll, Unnormalised()) for n > 1.
// At n = 1 taylorwin evaluates its formula, where Taylor returns [1], as every window here and in SciPy does.
//
// Gaussian takes the standard deviation in samples, same as SciPy.
// MATLAB's gausswin(n, alpha) is Gaussian(n, (n-1)/(2*alpha)).
//
// Tukey's alpha is unrestricted, same as SciPy and MATLAB:
// at or below 0 is Boxcar, and at or above 1 is Hann.
// Parameters with no meaning outside their range are checked where SciPy doesn't check them:
//   - std, p, sig, tau, sll and attenuation must be positive and finite
//   - nbar must be at least 1
//   - Kaiser's beta and the cosine sums' coefficients must be finite
//   - GeneralCosine needs a coefficient
//   - KaiserBesselDerived needs an even n, so n = 1 is an error rather than [1], and is symmetric only, same as SciPy
//
// # Algorithms
//
// Kaiser and KaiserBesselDerived use the modified Bessel function I0,
// ported from Cephes via SciPy's xsf, as THIRD_PARTY_NOTICES records.
// Chebwin is the real part of the FFT of its Chebyshev polynomial, through tukey/fft.
//
// # Conventions
//
// Every function is stateless:
// a new slice comes back, inputs are never modified, and all functions are safe for concurrent use.
// Invalid input returns an error wrapping ErrInvalidLength, ErrInvalidParameter or ErrInvalidOption,
// so check for one with errors.Is.
//
// A NaN or infinite parameter returns ErrInvalidParameter.
// The one exception is Tukey's alpha, where at or below 0 is Boxcar and at or above 1 is Hann, infinities included.
// For a flat window, use Boxcar rather than an infinite std or tau.
package window
