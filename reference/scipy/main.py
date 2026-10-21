"""
Reference vectors for tukey's conformance tests, generated with SciPy

Run with `just vectors-scipy`:
clears testdata/scipy first so renamed or removed cases get cleaned
then runs `uv run main.py <dir>` from reference/scipy

Requires pinned scipy/numpy in uv.lock: do not run against an unpinned environment

Writes <dir>/<category>/<name>.json, with one package per tukey GO package

Each JSON file records the scipy/numpy version
so a stale vector is very obvious in vc if run against different versions
"""

import argparse
from pathlib import Path

import dsp
import fft
import window
import write

# every generator in all generators
GENERATORS = dsp.GENERATORS + fft.GENERATORS + window.GENERATORS


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
