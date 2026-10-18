## The python version of the Survey Simulator

The python implementation builds an `ossssimlib` extension from the `F95`
sources (via f90wrap) and wraps it with Python classes in `ossssim`
(`OSSSSim`, model helpers, etc.).

See `pydoc ossssim` and the top-level [`examples/`](../examples/) directory.

Survey characterization formats and a mini example are documented in
[`../docs/`](../docs/). Model inputs and example characterizations are in
[`../data/`](../data/) (`Models/`, `Surveys/`); they are not shipped inside
this Python package.

### Installation

```bash
pip install .
```

(from the repository root). Requires gfortran and make. Pip’s isolated build
pulls f90wrap, NumPy, meson, and ninja from `[build-system]` in
`pyproject.toml`. The build runs `make -C F95 MODULE=ossssimlib` and installs
`_ossssimlib` plus the `ossssimlib` Python package alongside `ossssim`.

### Contents

- `ossssim` — Python API
- `ossssimlib` / `_ossssimlib` — f90wrap Fortran extension (built at install)
- `../tests` — unit tests
- `../examples` — usage examples
