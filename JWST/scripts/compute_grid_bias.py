#!/usr/bin/env python3
"""Grid-cell debiasing for JWST Sample A (Kavelaars et al. 2022 ac2c72 method)."""
from __future__ import annotations

import sys
from pathlib import Path

_SRC = Path(__file__).resolve().parents[2] / "src"
if str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))

from ossssim.grid_bias import A_STEP, H_STEP, JWST_SAMPLE_A, Q_STEP, SI_STEP
from ossssim.grid_bias_run import main

_HEADER = f"""# File: JWST-free-cla_m.detections-full
#
# Grid debiasing ac2c72; Eduardo et al. 2026 Sample A (20 objects)
# H_r from m_F150W2 + 1.0 - 5log10(r Δ) + 2.5log10(Bowell Φ), r=Δ=d_bary, G=-0.12
# Catalog i_free uses Ω=0 (node not in Sample A CSV)
#
# Grid size:
# h_step:  {H_STEP}
# a_step:  {A_STEP}
# q_step:  {Q_STEP}
# si_step: {SI_STEP}
#
"""


if __name__ == "__main__":
    main(
        JWST_SAMPLE_A,
        default_root=Path(__file__).resolve().parents[1],
        extra_header=_HEADER,
    )
