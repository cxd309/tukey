"""
Reference vectors for tukey's conformance tests, generated with SciPy

Run with `just vectors-scipy`:
clears testdata/scipy first so renamed or removed cases don't linger
then runs `uv run python -m vectors <dir>` from reference/scipy

Requires pinned scipy/numpy in uv.lock: do not run against an unpinned environment

Writes <dir>/<category>/<name>.json

Each JSON file records its own scipy/numpy version
so a stale vector is easy to spot even if this is re-run later against a newer pin
"""
