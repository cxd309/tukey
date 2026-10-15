"""Writing reference vector files to testdata/scipy"""

import json
from collections import Counter
from pathlib import Path

import numpy as np
import scipy as sc

# set from the command line by __main__ before any generator runs
output_root: Path | None = None
VERSIONS = {"scipy_version": sc.__version__, "numpy_version": np.__version__}

# vectors written per category this run, for the summary
written: Counter[str] = Counter()
# every file written this run, to catch two cases with the same name
seen: set[Path] = set()


def write_vector(
    category: str, name: str, description: str, params: dict, output: dict
) -> None:
    if output_root is None:
        raise RuntimeError("output_root not set: run via `just vectors`")
    out_dir = output_root / category
    out_dir.mkdir(parents=True, exist_ok=True)
    payload = {
        "description": description,
        **VERSIONS,
        "params": params,
        "output": output,
    }
    path = out_dir / f"{name}.json"
    if path in seen:
        raise ValueError(f"two cases named {category}/{name}")
    seen.add(path)
    path.write_text(json.dumps(payload, indent=2) + "\n")
    written[category] += 1


def complex_json(values) -> dict:
    """Complex values as {"re":[...], "im":[...]}, since JSON has no complex type"""
    values = np.asarray(values, dtype=complex)
    return {"re": values.real.tolist(), "im": values.imag.tolist()}
