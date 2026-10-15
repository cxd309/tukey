"""Write every reference vector: uv run python -m vectors <output dir>"""

import argparse
from pathlib import Path

from . import write
from .design import generate_butter, generate_zpk2sos
from .fft import generate_fft, generate_fft_options
from .filtering import (
    generate_filter,
    generate_filtfilt,
    generate_lfilter_state,
    generate_lfilter_zi,
    generate_sos_filtering,
)
from .freqz import generate_freqz

# every generator, each writes one or more categories
GENERATORS = [
    generate_butter,
    generate_zpk2sos,
    generate_filter,
    generate_lfilter_zi,
    generate_lfilter_state,
    generate_filtfilt,
    generate_sos_filtering,
    generate_freqz,
    generate_fft,
    generate_fft_options,
]


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Write tukey's SciPy reference vectors"
    )
    parser.add_argument(
        "output", type=Path, help="directory for <category>/<name>.json"
    )
    args = parser.parse_args()
    write.output_root = args.output

    for generate in GENERATORS:
        generate()

    for category, count in sorted(write.written.items()):
        print(f"{category}: {count}")
    print(f"{write.written.total()} vectors in {write.output_root}")


if __name__ == "__main__":
    main()
