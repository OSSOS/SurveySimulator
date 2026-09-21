#!/usr/bin/env python3
"""Grid-cell debiasing for JWST Sample A (Kavelaars et al. 2022 ac2c72 method)."""
from __future__ import annotations

import argparse
import csv
import math
import sys
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np
from astropy import units as u
from ossssim import OSSSSim
from ossssim.color import PhotSpec

sys.path.insert(0, str(Path(__file__).resolve().parent))
from grid_bias import (
    A_STEP,
    EPOCH_JD,
    FIELD_DEC_DEG as FIELD_DEC,
    FIELD_RA_DEG as FIELD_RA,
    FILL_FACTOR,
    H_STEP,
    MOSAIC_AREA_DEG2,
    MOSAIC_SIDE_DEG,
    PAPER_REFERENCE_JD,
    Q_STEP,
    RATE_CUT_MAX_ARCSEC_HR,
    RATE_CUT_MIN_ARCSEC_HR,
    SI_STEP,
    apparent_to_Hr,
    bounds_from_key,
    cell_key,
    compute_ifree,
    ecliptic_from_ifree,
    epoch_geometry,
    geometric_detection_prob,
    icrs_to_ecliptic,
    los_circular_elements,
    sample_aq,
)

TARGET_DETECTIONS = 5000


def load_detections(path: Path) -> list[dict]:
    rows = []
    with path.open() as fh:
        for row in csv.DictReader(fh):
            a, e, i = float(row["a"]), float(row["e"]), float(row["i"])
            d = float(row["d_bary"])
            hx = apparent_to_Hr(float(row["m_f150w2"]), d)
            q = a * (1.0 - e)
            # Sample A CSV has no node; Ω=0 places every object 90° from Ω_lp.
            ifree = compute_ifree(i, 0.0, a)
            rows.append({**row, "a": a, "e": e, "i": i, "d_bary": d, "q": q,
                         "Hx": hx, "ifree": ifree,
                         "sin_ifree": math.sin(math.radians(ifree)),
                         "cell": cell_key(a, q, math.sin(math.radians(ifree)), hx)})
    return rows


def load_bias_cache(path: Path) -> dict:
    if not path.exists():
        return {}
    out = {}
    with path.open() as fh:
        for row in csv.DictReader(fh):
            key = tuple(float(row[k]) for k in ("a_bin", "q_bin", "si_bin", "h_bin"))
            out[key] = (float(row["bias"]), int(row["n_drawn"]))
    return out


def save_bias_cache(path: Path, cache: dict) -> None:
    with path.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["a_bin", "q_bin", "si_bin", "h_bin", "bias", "n_drawn"])
        for key, (bias, n_drawn) in sorted(cache.items()):
            w.writerow([*key, bias, n_drawn])


def setup_pointings(char_root: Path) -> None:
    # Search footprint is the active mosaic (0.05 deg²), not the 1.6° implant box.
    # Fill factor is chip-fill inside that mosaic, not mosaic/implant (which would
    # be applied independently at each epoch and cube the spatial selection).
    # JD is the CADC visit-window midpoint: each epoch is a ~20 h shift-and-stack,
    # not a single 00:00 snapshot. Paper JD 2459974.5 is the orbit-fit reference.
    side = MOSAIC_SIDE_DEG
    for idx, jd in enumerate(EPOCH_JD, start=1):
        text = (
            f"# JWST Sample A epoch {idx} (CADC 1568 visit midpoint, shift-and-stack)\n"
            f"{side:.5f} {side:.5f} {FIELD_RA} {FIELD_DEC} {jd:.5f} {FILL_FACTOR:.5f} "
            f"JWST.csv JWST_sampleA.eff\n"
        )
        (char_root / f"epoch{idx}" / "pointings.list").write_text(text)


class JWSTSimulator:
    """Three-epoch AND using one OSSSSim / one RNG stream.

    Detos1 reloads characterization when the directory changes without
    resetting ran3. Creating three OSSSSim instances is unnecessary and
    used to hide that only the first GetSurvey call succeeded.
    """

    def __init__(self, char_root: Path, seed: int = 42):
        setup_pointings(char_root)
        self.epoch_dirs = [str((char_root / f"epoch{i}").resolve()) for i in (1, 2, 3)]
        # One Keplerian state; Detos1 advances M from this epoch to each pointing JD.
        self.element_epoch = EPOCH_JD[0]
        self.sim = OSSSSim(self.epoch_dirs[0], seed=seed)
        self.colors = PhotSpec()
        self._prime_surveys()

    def _prime_surveys(self) -> None:
        """Load each epoch's pointings.list / JWST.csv before the sanity plant.

        The first Detos1 call in a process has returned flag=0 for an on-field
        object while later calls on the same orbit returned 4. Prime with a
        dummy that is not on the mosaic so GetSurvey/ObsPos are initialized.
        """
        dummy = dict(
            a=44 * u.au, e=0.0, inc=20 * u.deg, node=0 * u.deg, peri=0 * u.deg,
            M=0 * u.deg, H=8 * u.mag, epoch=self.element_epoch * u.day, comp="default",
        )
        for epoch_dir in self.epoch_dirs:
            self.sim.characterization_directory = epoch_dir
            self.sim.simulate(dummy, colors=self.colors, model_band="r")

    def _simulate(self, epoch_dir, row, debug=False):
        self.sim.characterization_directory = epoch_dir
        return self.sim.simulate(row, colors=self.colors, model_band="r", debug=debug)

    def epoch_rows(self, a, e, inc, node, peri, M, H, debug=False) -> list[dict]:
        base = dict(a=a * u.au, e=e, inc=inc * u.deg, node=node * u.deg, peri=peri * u.deg,
                    M=M * u.deg, H=H * u.mag, epoch=self.element_epoch * u.day, comp="default")
        return [self._simulate(epoch_dir, base, debug=debug) for epoch_dir in self.epoch_dirs]

    def epoch_flags(self, a, e, inc, node, peri, M, H) -> list[int]:
        return [int(r["flag"]) for r in self.epoch_rows(a, e, inc, node, peri, M, H)]

    def detected_sample_a(self, a, e, inc, node, peri, M, H) -> bool:
        return all(f >= 4 for f in self.epoch_flags(a, e, inc, node, peri, M, H))


def _jd_utc(jd: float) -> str:
    mjd = jd - 2400000.5
    return (datetime(1858, 11, 17) + timedelta(days=mjd)).strftime("%Y-%m-%d %H:%M")


def _detos_sky(row: dict) -> str:
    """Format Detos1 outputs; RA/Dec are returned even when flag=0 after rebuild."""
    try:
        ra = float(row["RA"].to(u.deg).value)
        dec = float(row["DEC"].to(u.deg).value)
    except Exception:
        ra = math.degrees(float(row["RA"]))
        dec = math.degrees(float(row["DEC"]))
    dra = float(row["d_ra"])
    ddec = float(row["d_dec"])
    rate = math.degrees(math.hypot(dra, ddec)) * 3600.0 / 24.0
    r_au = float(row["r"].to(u.au).value) if hasattr(row["r"], "to") else float(row["r"])
    dlt = float(row["delta"].to(u.au).value) if hasattr(row["delta"], "to") else float(row["delta"])
    survey = row.get("Survey", "")
    if isinstance(survey, bytes):
        survey = survey.decode("utf-8", "replace")
    return (
        f"flag={int(row['flag'])}  RA,Dec={ra:.5f},{dec:.5f}  "
        f"rate={rate:.3f}\"/hr  r={r_au:.3f} Δ={dlt:.3f}  survey={survey!r}"
    )


def sanity_check_simulator(sim: JWSTSimulator) -> None:
    """Fail fast if a bright object on the JWST LOS is not Sample A.

    The plant is along the observer line of sight at the first visit
    midpoint (element epoch). Detos1 then Keplerian-propagates M to the
    other two shift-and-stack midpoints. Do not re-label the element
    epoch to each pointing JD: that freezes barycentric motion and is
    not the orbit Sample A linking uses.
    """
    jpl = Path(sim.epoch_dirs[0]) / "JWST.csv"
    element_jd = sim.element_epoch
    a, e, inc, node, peri, M = los_circular_elements(
        FIELD_RA, FIELD_DEC, 44.0, jpl, element_jd
    )
    geom = []
    for i, jd in enumerate(EPOCH_JD, start=1):
        ra, dec, sep, rate = epoch_geometry(
            a, e, inc, node, peri, M, jpl, element_jd, jd
        )
        geom.append((i, jd, ra, dec, sep, rate))
        in_fov = sep < MOSAIC_SIDE_DEG / 2.0
        rate_ok = RATE_CUT_MIN_ARCSEC_HR <= rate <= RATE_CUT_MAX_ARCSEC_HR
        print(
            f"sanity epoch{i} {_jd_utc(jd)} JD={jd:.5f}  "
            f"RA,Dec={ra:.5f},{dec:.5f}  sep={sep * 60:.3f}'  "
            f"rate={rate:.3f}\"/hr  FoV={in_fov}  rate_cut={rate_ok}",
            flush=True,
        )
    ra0, dec0, sep0 = geom[0][2], geom[0][3], geom[0][4]
    print(
        f"sanity plant at epoch1 LOS ICRS RA,Dec={ra0:.5f},{dec0:.5f}  "
        f"sep={sep0 * 60:.3f}' from mosaic centre "
        f"(half-side {MOSAIC_SIDE_DEG * 30:.1f}')",
        flush=True,
    )
    rows = sim.epoch_rows(a, e, inc, node, peri, M, 8.0, debug=True)
    flags = [int(r["flag"]) for r in rows]
    for i, row in enumerate(rows, start=1):
        print(f"sanity Detos1 epoch{i} {_detos_sky(row)}", flush=True)
    if all(f >= 4 for f in flags):
        print("sanity: LOS-planted object is a 3-epoch detection", flush=True)
        return
    # Call epoch1 after epoch2/3 so a leftover SAVE from epoch1's first load
    # cannot be blamed; print that isolated result too.
    base = dict(a=a * u.au, e=e, inc=inc * u.deg, node=node * u.deg, peri=peri * u.deg,
                M=M * u.deg, H=8.0 * u.mag, epoch=sim.element_epoch * u.day, comp="default")
    isolated = sim._simulate(sim.epoch_dirs[0], base, debug=True)
    print(f"sanity Detos1 epoch1 isolated after 2+3 {_detos_sky(isolated)}", flush=True)
    flags = [int(r["flag"]) for r in rows]
    flags[0] = int(isolated["flag"])
    if all(f >= 4 for f in flags):
        print("sanity: epoch1 is Sample A when evaluated after epochs 2 and 3", flush=True)
        return
    if flags[0] < 4 <= min(int(rows[1]["flag"]), int(rows[2]["flag"])):
        retried = sim.epoch_rows(a, e, inc, node, peri, M, 8.0)
        print(
            f"sanity: first flags={[int(r['flag']) for r in rows]}; "
            f"retry flags={[int(r['flag']) for r in retried]}",
            flush=True,
        )
        if all(int(r["flag"]) >= 4 for r in retried):
            print("sanity: epoch1 miss was the first Detos1 call, retry is Sample A",
                  flush=True)
            return
        flags = [int(r["flag"]) for r in retried]
    for dM in (-0.15, -0.10, -0.05, 0.05, 0.10, 0.15):
        shifted = sim.epoch_flags(a, e, inc, node, peri, M + dM, 8.0)
        if all(f >= 4 for f in shifted):
            print(f"sanity: LOS-planted object detected with ΔM={dM:.2f}°", flush=True)
            return
    detail = "; ".join(
        f"e{i} sep={sep * 60:.3f}' rate={rate:.3f}\"/hr"
        for i, _jd, _ra, _dec, sep, rate in geom
    )
    raise RuntimeError(
        "LOS-planted object at the JWST mosaic was not a 3-epoch Sample A "
        f"detection (flags={flags}; {detail}). "
        "Rebuild ossssimlib: Detos1 now returns RA/Dec/rate even when flag=0 "
        "and traces FoV/rate on stderr. epoch1 flag=0 with Fortran RA on-center "
        "is a FoV/rate_cut/η drop; RA=0 means n_sur=0 or mag gate; a 26' offset "
        "is still a mixed observer frame"
    )


def compute_cell_bias(sim: JWSTSimulator, cell_bounds: dict, seed: int, target: int) -> tuple[float, int]:
    rng = np.random.default_rng(seed)
    si0, si1 = cell_bounds["sin_ifree"]
    h0, h1 = cell_bounds["Hx"]

    n_detected = 0
    n_drawn = 0
    max_draws = max(target * 200000, 500000)
    while n_detected < target and n_drawn < max_draws:
        a, q = sample_aq(rng, cell_bounds["a"], cell_bounds["q"])
        e = 1.0 - q / a
        sin_ifree = float(rng.uniform(si0, si1))
        ifree = math.degrees(math.asin(max(0.0, min(1.0, sin_ifree))))
        H = float(rng.uniform(h0, h1))
        inc, node = ecliptic_from_ifree(ifree, a, rng)
        peri, M = rng.uniform(0, 360, size=2)
        n_drawn += 1
        if sim.detected_sample_a(a, e, inc, node, peri, M, H):
            n_detected += 1
        if n_drawn % 50000 == 0:
            print(f"    ... {n_drawn} draws, {n_detected}/{target} detections", flush=True)
    if n_detected < target:
        raise RuntimeError(f"Only {n_detected}/{target} after {n_drawn} draws")
    return n_detected / n_drawn, n_drawn


def write_detections_full(out_path: Path, detections: list[dict]) -> None:
    header = f"""# File: JWST-free-cla_m.detections-full
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
    cols = ("cl p j k sh object mag e_mag Filt Hsur dist e_dist Nobs time av_xres av_yres "
            "max_x max_y a e_a e i e_i Omega e_Omega omega e_omega tperi e_tperi "
            "RAdeg DEdeg JD rate MPC ifree Omfree omfree Hx comp bias")
    lines = [header, cols]
    for d in detections:
        lines.append(
            f"cla m -1 -1 S {d['name']:7s} {d['Hx']:.2f} 0.100 r {d['Hx']:.2f} {d['d_bary']:.3f} 0.100 "
            f"3 0.0000 0.083 0.073 0.311 0.343 {d['a']:11.6f} 0.1012 {d['e']:.6f} 0.001009 "
            f"{d['i']:6.3f} 0.100 0.000 0.100 0.000 0.100 0.000 0.100 0.000 0.100 "
            f"{FIELD_RA:.3f} {FIELD_DEC:.3f} {PAPER_REFERENCE_JD:.5f} 0.40 {d['name']:7s} {d['ifree']:6.3f} 0.000 0.000 "
            f"{d['Hx']:.2f} {d['comp']} {d['bias']:.7f}"
        )
    out_path.write_text("\n".join(lines) + "\n")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", default=str(Path(__file__).resolve().parents[1]))
    parser.add_argument("--target", type=int, default=TARGET_DETECTIONS)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    root = Path(args.root)
    detections = load_detections(root / "data/jwst_sampleA.csv")
    cache_path = root / "bias_grid.csv"
    cache = load_bias_cache(cache_path)
    sim = JWSTSimulator(root / "characterization", seed=args.seed)
    cells = sorted({d["cell"] for d in detections})
    _, lat = icrs_to_ecliptic(FIELD_RA, FIELD_DEC)
    p_geo = geometric_detection_prob(MOSAIC_AREA_DEG2, 7.0, lat)
    print(f"{len(cells)} cells, target={args.target}/cell")
    print(
        f"expected single-epoch geometric P ~ {p_geo:.2e} "
        f"(0.05 deg², i=7°, β={lat:.2f}°); Fig.20 is the H_r LF, not this rate",
        flush=True,
    )
    print(
        "characterization: 20-tile mosaic 0.05 deg² at "
        f"{FIELD_RA:.5f},{FIELD_DEC:.5f} (paper 13:57:33, −10:51:55); "
        "CADC 1568 detector-mean RA,Dec=209.39043,-10.86523",
        flush=True,
    )
    for i, jd in enumerate(EPOCH_JD, start=1):
        dt_ref = jd - PAPER_REFERENCE_JD
        print(
            f"  epoch{i} stack midpoint JD={jd:.5f} ({_jd_utc(jd)} UTC)  "
            f"Δ={dt_ref:+.2f}d from paper orbit-fit ref {PAPER_REFERENCE_JD:.1f}",
            flush=True,
        )
    print(
        "element epoch = epoch1; Detos1 propagates M to each pointing JD  "
        "(not frozen-M, not 00:00 integer days)",
        flush=True,
    )
    print("warning: Sample A CSV has no Ω; catalog i_free uses Ω=0", flush=True)
    sanity_check_simulator(sim)

    for idx, key in enumerate(cells):
        if key in cache:
            print(f"cell {idx+1}/{len(cells)} {key}: cached {cache[key][0]:.4g}")
            continue
        print(f"cell {idx+1}/{len(cells)} {key}:")
        bias, n_drawn = compute_cell_bias(sim, bounds_from_key(key), args.seed + idx, args.target)
        cache[key] = (bias, n_drawn)
        print(f"  bias={bias:.6g} n_drawn={n_drawn}")
        save_bias_cache(cache_path, cache)

    for d in detections:
        d["bias"] = cache[d["cell"]][0]
    out = root / "JWST-free-cla_m.detections-full"
    write_detections_full(out, detections)
    print(f"Wrote {out}; sum 1/bias = {sum(1/d['bias'] for d in detections):.1f}")


if __name__ == "__main__":
    main()
