#!/usr/bin/env python3
"""Grid-cell debiasing for JWST Sample A (Kavelaars et al. 2022 ac2c72 method)."""
from __future__ import annotations

import argparse
import csv
import math
import sys
from pathlib import Path

import numpy as np
from astropy import units as u
from ossssim import OSSSSim
from ossssim.color import PhotSpec

sys.path.insert(0, str(Path(__file__).resolve().parent))
from grid_bias import (
    A_STEP,
    FILL_FACTOR,
    H_STEP,
    MOSAIC_AREA_DEG2,
    MOSAIC_SIDE_DEG,
    Q_STEP,
    SI_STEP,
    apparent_to_Hr,
    apparent_radec_deg,
    bounds_from_key,
    cell_key,
    compute_ifree,
    ecliptic_from_ifree,
    geometric_detection_prob,
    icrs_to_ecliptic,
    los_circular_elements,
    parse_jpl_horizons_icrf,
    sample_aq,
    sky_separation_deg,
)

TARGET_DETECTIONS = 5000
EPOCH_JD = [2459969.5, 2459974.5, 2459978.5]
FIELD_RA = 209.3875
FIELD_DEC = -10.865278


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
    side = MOSAIC_SIDE_DEG
    for idx, jd in enumerate(EPOCH_JD, start=1):
        text = (
            f"# JWST Sample A epoch {idx}\n"
            f"{side:.5f} {side:.5f} {FIELD_RA} {FIELD_DEC} {jd} {FILL_FACTOR:.5f} "
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
        self.sim = OSSSSim(self.epoch_dirs[0], seed=seed)
        self.colors = PhotSpec()

    def epoch_flags(self, a, e, inc, node, peri, M, H) -> list[int]:
        base = dict(a=a * u.au, e=e, inc=inc * u.deg, node=node * u.deg, peri=peri * u.deg,
                    M=M * u.deg, H=H * u.mag, comp="default")
        flags = []
        for epoch_dir, jd in zip(self.epoch_dirs, EPOCH_JD):
            self.sim.characterization_directory = epoch_dir
            r = self.sim.simulate({**base, "epoch": jd * u.day}, colors=self.colors, model_band="r")
            flags.append(int(r["flag"]))
        return flags

    def detected_sample_a(self, a, e, inc, node, peri, M, H) -> bool:
        return all(f >= 4 for f in self.epoch_flags(a, e, inc, node, peri, M, H))


def sanity_check_simulator(sim: JWSTSimulator) -> None:
    """Fail fast if a bright object on the JWST LOS is not Sample A.

    The plant is along the observer line of sight, not the barycentric
    RA/Dec of the field (JWST parallax at 44 au is ~1°, larger than the mosaic).
    """
    jpl = Path(sim.epoch_dirs[0]) / "JWST.csv"
    a, e, inc, node, peri, M = los_circular_elements(
        FIELD_RA, FIELD_DEC, 44.0, jpl, EPOCH_JD[0]
    )
    obs = parse_jpl_horizons_icrf(jpl, EPOCH_JD[0])
    ra_pred, dec_pred = apparent_radec_deg(a, e, inc, node, peri, M, obs)
    sep = sky_separation_deg(ra_pred, dec_pred, FIELD_RA, FIELD_DEC)
    print(
        f"sanity plant ICRS RA,Dec={ra_pred:.5f},{dec_pred:.5f}  "
        f"sep={sep * 60:.3f}' from mosaic centre "
        f"(half-side {MOSAIC_SIDE_DEG * 30:.1f}')",
        flush=True,
    )
    flags = sim.epoch_flags(a, e, inc, node, peri, M, 8.0)
    if all(f >= 4 for f in flags):
        print("sanity: LOS-planted object is a 3-epoch detection", flush=True)
        return
    for dM in (-0.15, -0.10, -0.05, 0.05, 0.10, 0.15):
        shifted = sim.epoch_flags(a, e, inc, node, peri, M + dM, 8.0)
        if all(f >= 4 for f in shifted):
            print(f"sanity: LOS-planted object detected with ΔM={dM:.2f}°", flush=True)
            return
    raise RuntimeError(
        "LOS-planted object at the JWST mosaic was not a 3-epoch Sample A "
        f"detection (flags={flags}, predicted sep={sep * 60:.3f}'); "
        "Detos1 is not using the same observer frame as the plant "
        "(object ecliptic vs observatory ecliptic subtracted in ICRS)"
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
            f"{FIELD_RA:.3f} {FIELD_DEC:.3f} {EPOCH_JD[1]:.5f} 0.40 {d['name']:7s} {d['ifree']:6.3f} 0.000 0.000 "
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
