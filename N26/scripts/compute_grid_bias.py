#!/usr/bin/env python3
"""Grid-cell debiasing for Napier et al. 2026 N26 heliostack detections."""
from __future__ import annotations

import sys
from pathlib import Path

_SRC = Path(__file__).resolve().parents[2] / "src"
if str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))

from ossssim.grid_bias import A_STEP, H_STEP, N26_HELIOSTACK, Q_STEP, SI_STEP
from ossssim.grid_bias_run import main

_HEADER = f"""# File: N26-free-cla_m.detections-full
#
# Grid debiasing ac2c72; Napier et al. 2026 (PSJ 7, 117) heliostack of
# Bernstein 2004 GO-9433 ACS WFC F606W. Four objects (BF91, BH91, ABCD, WXYZ).
# r_AB = STMAG_F606W - 0.3; H_r from r_AB - 5log10(r Δ) + 2.5log10(Bowell Φ)
# a, e are Napier Table 2 a_min, e_min (CC-assumption lower bounds, 15-day arc)
# i_free is the midpoint of the Table 2 95% CI under the CC assumption
# Extra after MPC: ifree Omfree omfree (Laplace-free; Omfree=omfree=0), Hx, comp, bias
# Single 15-day stack (n_epochs=1), not 3-epoch AND
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
        N26_HELIOSTACK,
        default_root=Path(__file__).resolve().parents[1],
        extra_header=_HEADER,
    )
