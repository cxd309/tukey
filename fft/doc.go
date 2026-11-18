// Package fft computes discrete Fourier transforms of any length, in O(n log n),
// using only the standard library.
// It's checked against scipy.fft and MATLAB, which the tables in docs/accuracy.md measure.
//
//	y, err := fft.ForwardReal(x)
//	f, err := fft.RealFreq(len(x), 1/fs) // the frequency of each bin, in Hz
//
// # Functions
//
//	tukey          SciPy                    MATLAB
//	Forward        fft                      fft
//	Inverse        ifft                     ifft
//	ForwardReal    rfft                     fft, its first n/2 + 1 bins
//	InverseReal    irfft                    ifft(x, n, 'symmetric')
//	Freq           fftfreq
//	RealFreq       rfftfreq
//	Shift          fftshift                 fftshift
//	InverseShift   ifftshift                ifftshift
//	FastLength     next_fast_len(n, real=True)
//
// Options apply to every transform:
//   - NormBackward, the default, NormOrtho and NormForward choose where the 1/n scaling goes,
//     same as SciPy's norm
//   - N truncates or zero-pads the input first, same as SciPy's n and MATLAB's fft(x, n)
//
// # Conventions
//
// Every function is stateless:
// a new slice comes back, inputs are never modified, and all functions are safe for concurrent use.
// Invalid input returns an error wrapping ErrInvalidLength, ErrSignalTooShort,
// ErrInvalidParameter or ErrInvalidOption, so check for one with errors.Is.
//
// A NaN or infinite sample is never an error:
// A transform mixes every sample into every bin, so it makes every bin NaN.
// A NaN or infinite parameter, such as Freq's d, returns ErrInvalidParameter.
//
// # Where tukey differs
//
//   - the names read after the package, so fft.Forward rather than fft.FFT
//   - InverseReal needs its length, as n/2 + 1 bins come from both an even and an odd length
//   - Freq lists the upper half of the bins as negative frequencies, same as SciPy,
//     where MATLAB's usual (0:n-1)*fs/n lists them above Nyquist
//   - FastLength allows factors of 2, 3 and 5 only, where SciPy's next_fast_len also allows 7 and 11
//   - an empty input without N is ErrSignalTooShort, same as SciPy
//
// # Performance
//
// Every length is O(n log n).
// Powers of two use iterative radix-2 Cooley-Tukey,
// other lengths recursive mixed radix with dedicated butterflies for 2, 3, 4 and 5,
// and primes above 13 Bluestein's algorithm.
// Lengths with only 2, 3 and 5 as factors are the fastest, so zero-pad to FastLength when any length will do.
// A real transform of even length costs about half a complex one.
// Twiddle factors are reduced to an eighth of a turn exactly, so quarter and half turns are exact.
package fft
