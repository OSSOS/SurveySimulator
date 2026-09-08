## The python version of the Survey Simulator

The python implementation builds an `ossssimlib` extension from the `F95`
sources (via f90wrap) and wraps it with Python classes in `ossssim`
(`OSSSSim`, model helpers, etc.).

See `pydoc ossssim` and the top-level [`examples/`](../examples/) directory.

Survey characterization formats and a mini example are documented in
[`../docs/`](../docs/). Full survey characterizations and model files are
distributed separately as [`SurveySimulator-Data`](../../SurveySimulator-Data/)
(DOI / web download; not shipped inside this Python package). In-repo test
fixtures live under `../F95/tests/Surveys/` and `../F95/tests/Models/`.

### Installation

```bash
pip install .
```

(from the repository root). Requires gfortran, make, and f90wrap. The build
runs `make -C F95 MODULE=ossssimlib` and installs `_ossssimlib` plus the
`ossssimlib` Python package alongside `ossssim`.

### Contents

- `ossssim` — Python API
- `ossssimlib` / `_ossssimlib` — f90wrap Fortran extension (built at install)
- `../tests` — unit tests
- `../examples` — usage examples
