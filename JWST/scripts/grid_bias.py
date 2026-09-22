"""JWST Sample A grid-cell helpers.

Loads ``src/ossssim/grid_bias.py`` without importing the ``ossssim`` package
(that would pull ossssimlib). Module-level field defaults are this survey.
"""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

_PATH = Path(__file__).resolve().parents[2] / "src" / "ossssim" / "grid_bias.py"
_NAME = "ossssim_grid_bias_impl"
if _NAME not in sys.modules:
    spec = importlib.util.spec_from_file_location(_NAME, _PATH)
    _mod = importlib.util.module_from_spec(spec)
    sys.modules[_NAME] = _mod
    spec.loader.exec_module(_mod)
else:
    _mod = sys.modules[_NAME]

globals().update({k: v for k, v in vars(_mod).items() if not k.startswith("__")})
SURVEY = JWST_SAMPLE_A  # noqa: F405
