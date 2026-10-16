// Package fft computes discrete Fourier transforms of any length
// matches scipy.fft, standard library only
//
// # Transforms
//
// # Forward and Inverse transform complex signals
//
// ForwardReal and InverseReal transform real ones,
// keeping only the n/2 + 1 non-negative frequency bins,
// since a real signal's spectrum mirrors itself
//
// Freq and RealFreq give each bin's frequency,
// Shift and InverseShift move zero frequency to the centre and back, for plotting
//
// Options apply to every transform:
// Norm chooses where the 1/n scaling goes
// NormBackward is the default but also support NormOrtho and NormForward
//
// N truncates or zero-pads the input first, as MATLAB's fft(x, n).
//
// # Compatibility
//
// Every function matches its scipy.fft counterpart within a measured tolerance,
// checked against reference vectors at 23 lengths covering every algorithm path:
//
//	Function       SciPy equivalent         Worst measured   Tested to
//	Forward        fft                      9.3e-16          1e-14
//	Inverse        ifft                     8.5e-16          1e-14
//	ForwardReal    rfft                     5.8e-16          1e-14
//	InverseReal    irfft                    5.0e-16          1e-14
//	Norm, N        norm=, n=                6.4e-16          1e-14
//	Freq           fftfreq, rfftfreq        exact            exact
//	Shift          fftshift, ifftshift      exact            exact
//
// Errors are relative to the largest expected value
// 1e-16 is about 1ULP of a float64
//
// Every transform is also checked against MATLAB R2025a, whose fft is FFTW,
// an independent implementation, on the same inputs:
//
//	Function       MATLAB equivalent        Worst measured   Tested to
//	Forward        fft                      9.5e-16          1e-14
//	Inverse        ifft                     8.3e-16          1e-14
//	ForwardReal    fft, first n/2+1 bins    6.7e-16          1e-14
//	InverseReal    ifft(x, n, 'symmetric')  5.3e-16          1e-14
//	N              fft(x, n), ifft(x, n)    7.2e-16          1e-14
//	Shift          fftshift, ifftshift      exact            exact
//
// MATLAB has no normalisation option or fftfreq, so Norm and Freq are checked
// against SciPy only
//
// # Algorithms
//
// Lengths that are powers of two use iterative radix-2 Cooley–Tukey
// Other lengths use recursive mixed radix
//   - dedicated butterflies factors of 2,3,4,5
//   - direct sum for primes up to 13
//   - Bluestein's algorithm for larger primes
//
// every length is O(n log n)
//
// Real transforms of even length pack pairs of samples into one complex value
// and transform half the length
//
// Twiddle factors are computed by reducing each angle to an eighth of a turn with
// exact integer arithmetic, so quarter and half turns are exact.
//
// # Conventions
//
// Every function is stateless: a new slice comes back, inputs are never
// modified, and all functions are safe for concurrent use.
// Invalid input returns an error wrapping ErrInvalidLength or ErrInvalidOption,
// can easily be checked with errors.Is
package fft
