// Package window computes window functions, matching scipy.signal.windows,
// using only the standard library.
// It has every SciPy window except dpss,
// and is checked against SciPy and MATLAB, which the README's compatibility tables measure.
//
//	w, err := window.Hann(1024, window.Periodic()) // for spectral analysis
//	h, err := window.Hamming(101)                  // for filter design
//
// Each window returns n samples:
// symmetric by default, for filter design,
// or periodic, DFT-even, with Periodic, for spectral analysis.
//
// # Windows
//
//	tukey                SciPy                    MATLAB
//	Hann                 hann                     hann
//	Hamming              hamming                  hamming
//	Blackman             blackman                 blackman
//	BlackmanHarris       blackmanharris           blackmanharris
//	Nuttall              nuttall                  nuttallwin
//	FlatTop              flattop                  flattopwin
//	GeneralHamming       general_hamming
//	GeneralCosine        general_cosine
//	Boxcar               boxcar                   rectwin
//	Triang               triang                   triang
//	Bartlett             bartlett                 bartlett
//	Parzen               parzen                   parzenwin
//	Bohman               bohman                   bohmanwin
//	Barthann             barthann                 barthannwin
//	Tukey                tukey                    tukeywin
//	Cosine               cosine
//	Lanczos              lanczos
//	Gaussian             gaussian                 gausswin
//	GeneralGaussian      general_gaussian
//	Exponential          exponential
//	Kaiser               kaiser                   kaiser
//	KaiserBesselDerived  kaiser_bessel_derived
//	Taylor               taylor                   taylorwin
//	Chebwin              chebwin                  chebwin
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
//
// # Options
//
// One Option type serves every window, and a window given an option that doesn't apply to it returns ErrInvalidOption:
//   - Periodic gives the symmetric window of n+1 points without its last point, for every window but KaiserBesselDerived.
//     It's SciPy's sym=False and MATLAB's 'periodic'.
//   - Center moves Exponential's peak, with Periodic only.
//   - Unnormalised leaves Taylor's peak where its formula puts it, rather than at 1.
//     It's SciPy's norm=False, and MATLAB's taylorwin.
//
// # Where tukey differs
//
//   - symmetric windows are computed for their first half and mirrored, so they're exactly symmetric,
//     where SciPy's are symmetric to within a few ULP
//   - MATLAB's rectwin, triang, bartlett, parzenwin, bohmanwin, barthannwin, tukeywin, gausswin, kaiser,
//     taylorwin and chebwin have no periodic form
//   - MATLAB's hanning is Hann without its zero ends, and Hann's doc shows how to build it
//   - Gaussian takes the standard deviation in samples, same as SciPy,
//     and MATLAB's gausswin(n, alpha) is Gaussian(n, (n-1)/(2*alpha))
//   - Taylor takes its sidelobe level as positive dB, same as SciPy,
//     and MATLAB's taylorwin(n, nbar, -sll) is Taylor(n, nbar, sll, Unnormalised()) for n > 1
//   - every window is [1] at n = 1, same as SciPy, where MATLAB's taylorwin evaluates its formula
//   - Tukey's alpha is unrestricted, same as SciPy and MATLAB, so at or below 0 is Boxcar, and at or above 1 is Hann
//   - parameters with no meaning outside their range are checked where SciPy doesn't check them:
//     std, p, sig, tau, sll and attenuation must be positive and finite, nbar at least 1,
//     Kaiser's beta and the cosine sums' coefficients finite,
//     and KaiserBesselDerived's n even, so n = 1 is an error rather than [1]
//   - Chebwin below about 45 dB isn't suitable for spectral analysis, as SciPy warns
package window
