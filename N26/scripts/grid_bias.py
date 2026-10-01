"""N26 heliostack grid-cell helpers (Napier et al. 2026).

Loads the shared ossssim.grid_bias implementation without importing the
ossssim package (that would pull ossssimlib), then binds N26_HELIOSTACK as
the default survey.
"""
from __future__ import annotations

import importlib.util
import sys
from functools import wraps
from pathlib import Path

_PATH = Path(__file__).resolve().parents[2] / "src" / "ossssim" / "grid_bias.py"
_NAME = "ossssim_grid_bias_impl"
if _NAME not in sys.modules:
    spec = importlib.util.spec_from_file_location(_NAME, _PATH)
    _gb = importlib.util.module_from_spec(spec)
    sys.modules[_NAME] = _gb
    spec.loader.exec_module(_gb)
else:
    _gb = sys.modules[_NAME]

globals().update({k: v for k, v in vars(_gb).items() if not k.startswith("__")})

SURVEY = _gb.N26_HELIOSTACK
FIELD_RA_DEG = SURVEY.field_ra_deg
FIELD_DEC_DEG = SURVEY.field_dec_deg
MOSAIC_AREA_DEG2 = SURVEY.mosaic_area_deg2
MOSAIC_SIDE_DEG = SURVEY.mosaic_side_deg
MOSAIC_WIDTH_DEG = SURVEY.mosaic_width_deg
MOSAIC_HEIGHT_DEG = SURVEY.mosaic_height_deg
FILL_FACTOR = SURVEY.fill_factor
EPOCH_JD = SURVEY.epoch_jd
PAPER_REFERENCE_JD = SURVEY.paper_reference_jd
RATE_CUT_MIN_ARCSEC_HR = SURVEY.rate_cut_min_arcsec_hr
RATE_CUT_MAX_ARCSEC_HR = SURVEY.rate_cut_max_arcsec_hr
H_COLOR_OFFSET = SURVEY.mag_color_offset


def _with_n26(fn):
    @wraps(fn)
    def wrapped(*args, survey=SURVEY, **kwargs):
        return fn(*args, survey=survey, **kwargs)
    return wrapped


apparent_to_Hr = _with_n26(_gb.apparent_to_Hr)
geometric_prob_for_aimed = _with_n26(_gb.geometric_prob_for_aimed)
aimed_at_field = _with_n26(_gb.aimed_at_field)
epoch_geometry = _with_n26(_gb.epoch_geometry)
render_pointings_text = _with_n26(_gb.render_pointings_text)
setup_pointings = _with_n26(_gb.setup_pointings)
sample_mosaic_icrs = _with_n26(_gb.sample_mosaic_icrs)
sample_aimed_elements = _with_n26(_gb.sample_aimed_elements)
write_bias_check_plots = _with_n26(_gb.write_bias_check_plots)
load_detections = _with_n26(_gb.load_detections)
write_detections_full = _with_n26(_gb.write_detections_full)
