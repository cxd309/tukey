"""Reference vectors for tukey/dsp: filter design, filtering and frequency response"""

from dsp.design import generate_butter, generate_zpk2sos
from dsp.filtering import (
    generate_filter,
    generate_filtfilt,
    generate_lfilter_state,
    generate_lfilter_zi,
    generate_sos_filtering,
)
from dsp.freqz import generate_freqz

# every dsp generator; each writes one or more categories
GENERATORS = [
    generate_butter,
    generate_zpk2sos,
    generate_filter,
    generate_lfilter_zi,
    generate_lfilter_state,
    generate_filtfilt,
    generate_sos_filtering,
    generate_freqz,
]
