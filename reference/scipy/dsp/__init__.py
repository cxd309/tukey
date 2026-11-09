"""Reference vectors for tukey/dsp: filter design, filtering and frequency response"""

from dsp.alignment import generate_findsignal
from dsp.conditioning import (
    generate_detrend,
    generate_fillmissing,
    generate_interp1,
    generate_rms,
)
from dsp.convolution import generate_convolve, generate_correlate
from dsp.design import (
    generate_butter,
    generate_cheby1,
    generate_cheby2,
    generate_zpk2sos,
)
from dsp.envelope import generate_envelope, generate_hilbert
from dsp.filtering import (
    generate_filter,
    generate_filtfilt,
    generate_lfilter_state,
    generate_lfilter_zi,
    generate_sos_filtering,
)
from dsp.freqz import generate_freqz
from dsp.peaks import (
    generate_find_peaks,
    generate_peak_prominences,
    generate_peak_widths,
)
from dsp.reductions import (
    generate_cumsum,
    generate_cumulative_trapezoid,
    generate_diff,
    generate_gradient,
    generate_mean,
    generate_percentile,
    generate_sum,
    generate_trapezoid,
)
from dsp.spectral import (
    generate_periodogram,
    generate_spectrogram,
    generate_stft,
    generate_stft_frequencies,
    generate_welch,
)

# every dsp generator; each writes one or more categories
GENERATORS = [
    generate_butter,
    generate_cheby1,
    generate_cheby2,
    generate_zpk2sos,
    generate_filter,
    generate_lfilter_zi,
    generate_lfilter_state,
    generate_filtfilt,
    generate_sos_filtering,
    generate_freqz,
    generate_detrend,
    generate_interp1,
    generate_fillmissing,
    generate_rms,
    generate_find_peaks,
    generate_peak_prominences,
    generate_peak_widths,
    generate_hilbert,
    generate_envelope,
    generate_findsignal,
    generate_stft,
    generate_stft_frequencies,
    generate_spectrogram,
    generate_welch,
    generate_periodogram,
    generate_sum,
    generate_mean,
    generate_diff,
    generate_gradient,
    generate_cumsum,
    generate_trapezoid,
    generate_cumulative_trapezoid,
    generate_percentile,
    generate_convolve,
    generate_correlate,
]
