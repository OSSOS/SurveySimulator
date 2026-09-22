"""Grid-cell Horvitz–Thompson runner shared by JWST Sample A and N26."""
from __future__ import annotations

import argparse
import csv
import math
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np
from astropy import units as u

from ossssim import OSSSSim
from ossssim.color import PhotSpec
from ossssim.grid_bias import (
    A_STEP,
    H_STEP,
    Q_STEP,
    SI_STEP,
    GridSurvey,
    JWST_SAMPLE_A,
    aimed_detection_bias,
    as_check_arrays,
    bounds_from_key,
    check_plot_tag,
    empty_check_samples,
    epoch_geometry,
    geometric_detection_prob,
    geometric_prob_for_aimed,
    icrs_to_ecliptic,
    load_detections,
    los_circular_elements,
    parse_jpl_horizons_icrf,
    record_check_sample,
    sample_aimed_elements,
    sample_aq,
    setup_pointings,
    stack_check_samples,
    write_bias_check_plots,
)

TARGET_DETECTIONS = 5000


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


class GridBiasSimulator:
    """One OSSSSim / one RNG stream over the survey's characterization epochs.

    Detos1 reloads characterization when the directory changes without
    resetting ran3. Creating one OSSSSim per epoch is unnecessary and
    used to hide that only the first GetSurvey call succeeded.

    JWST Sample A is 3-epoch AND (flag≥4 at every visit-stack). N26 is a
    single 15-day heliostack, so ``n_epochs=1``.
    """

    def __init__(self, survey: GridSurvey, char_root: Path, seed: int = 42):
        self.survey = survey
        setup_pointings(char_root, survey=survey)
        char_root = Path(char_root)
        if survey.epoch_layout == "flat":
            self.epoch_dirs = [str(char_root.resolve())]
        else:
            self.epoch_dirs = [
                str((char_root / f"epoch{i}").resolve())
                for i in range(1, survey.n_epochs + 1)
            ]
        self.element_epoch = survey.epoch_jd[0]
        self.sim = OSSSSim(self.epoch_dirs[0], seed=seed)
        self.colors = PhotSpec()
        self._prime_surveys()

    def _prime_surveys(self) -> None:
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

    def detected(self, a, e, inc, node, peri, M, H) -> bool:
        return all(f >= 4 for f in self.epoch_flags(a, e, inc, node, peri, M, H))


class JWSTSimulator(GridBiasSimulator):
    """JWST Sample A 3-epoch AND. ``char_root`` is the first positional arg."""

    def __init__(self, char_root: Path, seed: int = 42, survey: GridSurvey = JWST_SAMPLE_A):
        super().__init__(survey, char_root, seed)


def _jd_utc(jd: float) -> str:
    mjd = jd - 2400000.5
    return (datetime(1858, 11, 17) + timedelta(days=mjd)).strftime("%Y-%m-%d %H:%M")


def _row_radec(row: dict) -> tuple[float, float]:
    """Detos1 ICRS RA/Dec in degrees (available even when flag=0 after rebuild)."""
    try:
        ra = float(row["RA"].to(u.deg).value)
        dec = float(row["DEC"].to(u.deg).value)
    except Exception:
        ra = math.degrees(float(row["RA"]))
        dec = math.degrees(float(row["DEC"]))
    return ra, dec


def _detos_sky(row: dict) -> str:
    ra, dec = _row_radec(row)
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


def sanity_check_simulator(sim: GridBiasSimulator) -> None:
    """Fail fast if a bright object on the survey LOS is not a detection."""
    survey = sim.survey
    jpl = Path(sim.epoch_dirs[0]) / survey.observer_csv
    element_jd = sim.element_epoch
    a, e, inc, node, peri, M = los_circular_elements(
        survey.field_ra_deg, survey.field_dec_deg, 44.0, jpl, element_jd
    )
    geom = []
    half_w = 0.5 * survey.mosaic_width_deg
    half_h = 0.5 * survey.mosaic_height_deg
    for i, jd in enumerate(survey.epoch_jd, start=1):
        ra, dec, sep, rate = epoch_geometry(
            a, e, inc, node, peri, M, jpl, element_jd, jd, survey=survey
        )
        geom.append((i, jd, ra, dec, sep, rate))
        in_fov = (abs(ra - survey.field_ra_deg) <= half_w + 1e-3
                  and abs(dec - survey.field_dec_deg) <= half_h + 1e-3)
        rate_ok = (survey.rate_cut_min_arcsec_hr <= rate <= survey.rate_cut_max_arcsec_hr)
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
        f"({survey.mosaic_width_deg * 60:.1f}' × {survey.mosaic_height_deg * 60:.1f}')",
        flush=True,
    )
    rows = sim.epoch_rows(a, e, inc, node, peri, M, 8.0, debug=True)
    flags = [int(r["flag"]) for r in rows]
    for i, row in enumerate(rows, start=1):
        print(f"sanity Detos1 epoch{i} {_detos_sky(row)}", flush=True)
    if all(f >= 4 for f in flags):
        print(f"sanity: LOS-planted object is a {survey.n_epochs}-epoch detection",
              flush=True)
        return
    base = dict(a=a * u.au, e=e, inc=inc * u.deg, node=node * u.deg, peri=peri * u.deg,
                M=M * u.deg, H=8.0 * u.mag, epoch=sim.element_epoch * u.day, comp="default")
    isolated = sim._simulate(sim.epoch_dirs[0], base, debug=True)
    print(f"sanity Detos1 epoch1 isolated {_detos_sky(isolated)}", flush=True)
    flags = [int(r["flag"]) for r in rows]
    flags[0] = int(isolated["flag"])
    if all(f >= 4 for f in flags):
        print("sanity: epoch1 is a detection when evaluated after later epochs",
              flush=True)
        return
    if flags[0] < 4 and survey.n_epochs > 1 and min(flags[1:]) >= 4:
        retried = sim.epoch_rows(a, e, inc, node, peri, M, 8.0)
        print(
            f"sanity: first flags={[int(r['flag']) for r in rows]}; "
            f"retry flags={[int(r['flag']) for r in retried]}",
            flush=True,
        )
        if all(int(r["flag"]) >= 4 for r in retried):
            print("sanity: epoch1 miss was the first Detos1 call, retry is a detection",
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
        f"LOS-planted object at the {survey.name} mosaic was not a "
        f"{survey.n_epochs}-epoch detection (flags={flags}; {detail})."
    )


def compute_cell_bias(sim: GridBiasSimulator, cell_bounds: dict, seed: int, target: int,
                      plot_dir: Path | None = None, plot_tags: list | None = None
                      ) -> tuple[float, int, dict, dict]:
    """P(detect | cell) by FoV-aimed (Ω, ω, M) times single-epoch P_geom."""
    survey = sim.survey
    rng = np.random.default_rng(seed)
    si0, si1 = cell_bounds["sin_ifree"]
    h0, h1 = cell_bounds["Hx"]
    jpl = Path(sim.epoch_dirs[0]) / survey.observer_csv
    obs = parse_jpl_horizons_icrf(jpl, sim.element_epoch)
    sampled = empty_check_samples()
    detected = empty_check_samples()

    n_detected = 0
    n_aimed = 0
    n_fail = 0
    geom_weight_sum = 0.0
    max_tries = max(target * 1000, 10000)
    while n_detected < target and (n_aimed + n_fail) < max_tries:
        a, q = sample_aq(rng, cell_bounds["a"], cell_bounds["q"])
        e = 1.0 - q / a
        sin_ifree = float(rng.uniform(si0, si1))
        ifree = math.degrees(math.asin(max(0.0, min(1.0, sin_ifree))))
        H = float(rng.uniform(h0, h1))
        el = sample_aimed_elements(a, e, ifree, obs, rng, survey=survey)
        if el is None:
            n_fail += 1
            continue
        inc, node, peri, M = el
        n_aimed += 1
        p_geom = geometric_prob_for_aimed(a, e, inc, node, peri, M, survey=survey)
        rows = sim.epoch_rows(a, e, inc, node, peri, M, H)
        flags = [int(r["flag"]) for r in rows]
        ra, dec = _row_radec(rows[0])
        record_check_sample(sampled, ra, dec, a, e, inc, node, peri, M)
        if all(f >= 4 for f in flags):
            n_detected += 1
            geom_weight_sum += p_geom
            record_check_sample(detected, ra, dec, a, e, inc, node, peri, M)
        if n_aimed % 500 == 0:
            bias_so_far = aimed_detection_bias(n_aimed, geom_weight_sum)
            print(
                f"    ... {n_aimed} aimed ({n_fail} invert-fail), "
                f"{n_detected}/{target} detections  "
                f"P(det|FoV)={n_detected / n_aimed:.3g}  bias~{bias_so_far:.3g}",
                flush=True,
            )
    sampled_arr = as_check_arrays(sampled)
    detected_arr = as_check_arrays(detected)
    if plot_dir is not None and sampled_arr["ra"].size:
        for tag in plot_tags or []:
            for path in write_bias_check_plots(
                    plot_dir, sampled_arr, detected_arr, check_plot_tag(tag),
                    survey=survey,
            ):
                print(f"    wrote {path}", flush=True)
    if n_aimed == 0:
        return 0.0, 0, sampled_arr, detected_arr
    if n_detected < target:
        raise RuntimeError(
            f"Only {n_detected}/{target} after {n_aimed} aimed plants "
            f"({n_fail} invert-fail)"
        )
    return aimed_detection_bias(n_aimed, geom_weight_sum), n_aimed, sampled_arr, detected_arr


def write_detections_full(out_path: Path, detections: list[dict], survey: GridSurvey,
                          header_lines: str | None = None) -> None:
    if header_lines is None:
        header_lines = (
            f"# File: {survey.detections_full_name}\n"
            f"#\n"
            f"# Grid debiasing ac2c72; {survey.name}\n"
            f"# H_r from {survey.mag_column}{survey.mag_color_offset:+.1f} "
            f"- 5log10(r Δ) + 2.5log10(Bowell Φ), r=Δ=d_bary, G=-0.12\n"
            f"#\n"
            f"# Grid size:\n"
            f"# h_step:  {H_STEP}\n"
            f"# a_step:  {A_STEP}\n"
            f"# q_step:  {Q_STEP}\n"
            f"# si_step: {SI_STEP}\n"
            f"#\n"
        )
    cols = ("cl p j k sh object mag e_mag Filt Hsur dist e_dist Nobs time av_xres av_yres "
            "max_x max_y a e_a e i e_i Omega e_Omega omega e_omega tperi e_tperi "
            "RAdeg DEdeg JD rate MPC ifree Omfree omfree Hx comp bias")
    ref_jd = survey.paper_reference_jd or survey.epoch_jd[0]
    nobs = survey.n_epochs
    lines = [header_lines, cols]
    for d in detections:
        name = str(d["name"])
        lines.append(
            f"cla m -1 -1 S {name:7s} {d['Hx']:.2f} 0.100 r {d['Hx']:.2f} {d['d_bary']:.3f} 0.100 "
            f"{nobs} 0.0000 0.083 0.073 0.311 0.343 {d['a']:11.6f} 0.1012 {d['e']:.6f} 0.001009 "
            f"{d['i']:6.3f} 0.100 0.000 0.100 0.000 0.100 0.000 0.100 0.000 0.100 "
            f"{survey.field_ra_deg:.3f} {survey.field_dec_deg:.3f} {ref_jd:.5f} 0.40 {name:7s} {d['ifree']:6.3f} 0.000 0.000 "
            f"{d['Hx']:.2f} {d['comp']} {d['bias']:.7f}"
        )
    out_path.write_text("\n".join(lines) + "\n")


def run_grid_bias(survey: GridSurvey, root: Path, target: int = TARGET_DETECTIONS,
                  seed: int = 42, check_plots_dir: Path | None = None,
                  no_check_plots: bool = False,
                  extra_header: str | None = None) -> Path:
    detections = load_detections(root / survey.detections_relpath, survey)
    cache_path = root / "bias_grid.csv"
    cache = load_bias_cache(cache_path)
    sim = GridBiasSimulator(survey, root / "characterization", seed=seed)
    cells = sorted({d["cell"] for d in detections})
    _, lat = icrs_to_ecliptic(survey.field_ra_deg, survey.field_dec_deg)
    p_geo = geometric_detection_prob(survey.mosaic_area_deg2, 7.0, lat)
    print(f"{len(cells)} cells, target={target}/cell")
    print(
        f"expected single-epoch geometric P ~ {p_geo:.2e} "
        f"({survey.mosaic_area_deg2:.4f} deg², i=7°, β={lat:.2f}°)",
        flush=True,
    )
    print(
        "sampling: draw a/e/i_free/H in the cell, sample r on [q, Q], invert "
        "(Ω, ω, M) onto the ICRS mosaic (LOS rotated to ecliptic); HT bias is "
        "P(detect | FoV) × P_geom",
        flush=True,
    )
    print(
        f"characterization: {survey.mosaic_width_deg:.5f}×{survey.mosaic_height_deg:.5f} deg "
        f"({survey.mosaic_area_deg2:.4f} deg²) at "
        f"{survey.field_ra_deg:.5f},{survey.field_dec_deg:.5f}",
        flush=True,
    )
    ref_jd = survey.paper_reference_jd
    for i, jd in enumerate(survey.epoch_jd, start=1):
        extra = ""
        if ref_jd is not None:
            extra = f"  Δ={jd - ref_jd:+.2f}d from paper ref {ref_jd:.1f}"
        print(
            f"  epoch{i} stack midpoint JD={jd:.5f} ({_jd_utc(jd)} UTC){extra}",
            flush=True,
        )
    print(
        "element epoch = epoch1; Detos1 propagates M to each pointing JD",
        flush=True,
    )
    sanity_check_simulator(sim)

    plot_dir = None if no_check_plots else Path(
        check_plots_dir or (root / "check_plots")
    )
    if plot_dir is not None:
        plot_dir.mkdir(parents=True, exist_ok=True)
    print(f"check plots → {plot_dir} (one RA/Dec + elements pair per object)", flush=True)
    run_sampled = []
    run_detected = []

    for idx, key in enumerate(cells):
        if key in cache:
            print(f"cell {idx+1}/{len(cells)} {key}: cached {cache[key][0]:.4g}")
            continue
        print(f"cell {idx+1}/{len(cells)} {key}:")
        members = [d["name"] for d in detections if d["cell"] == key]
        print(f"  objects: {', '.join(str(n) for n in members)}", flush=True)
        bias, n_drawn, sampled, detected = compute_cell_bias(
            sim, bounds_from_key(key), seed + idx, target,
            plot_dir=plot_dir, plot_tags=members,
        )
        cache[key] = (bias, n_drawn)
        print(f"  bias={bias:.6g} n_drawn={n_drawn}")
        save_bias_cache(cache_path, cache)
        run_sampled.append(sampled)
        run_detected.append(detected)

    if plot_dir is not None and run_sampled:
        all_s = stack_check_samples(run_sampled)
        all_d = stack_check_samples(run_detected)
        if all_s["ra"].size:
            for path in write_bias_check_plots(plot_dir, all_s, all_d, "all",
                                               survey=survey):
                print(f"wrote {path} (all objects combined)", flush=True)

    for d in detections:
        d["bias"] = cache[d["cell"]][0]
    out = root / survey.detections_full_name
    write_detections_full(out, detections, survey, header_lines=extra_header)
    print(f"Wrote {out}; sum 1/bias = {sum(1/d['bias'] for d in detections):.1f}")
    return out


def build_arg_parser(survey: GridSurvey, default_root: Path) -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=f"Grid-cell debiasing for {survey.name} (Kavelaars et al. 2022 ac2c72)."
    )
    parser.add_argument("--root", default=str(default_root))
    parser.add_argument("--target", type=int, default=TARGET_DETECTIONS)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument(
        "--check-plots-dir", default=None,
        help="Directory for sampled-vs-detected RA/Dec and element check plots "
             "(default: <root>/check_plots)",
    )
    parser.add_argument(
        "--no-check-plots", action="store_true",
        help="Skip writing sampled-vs-detected check plots",
    )
    return parser


def main(survey: GridSurvey | None = None, default_root: Path | None = None,
         extra_header: str | None = None) -> None:
    survey = survey or JWST_SAMPLE_A
    default_root = default_root or Path.cwd()
    args = build_arg_parser(survey, default_root).parse_args()
    run_grid_bias(
        survey, Path(args.root), target=args.target, seed=args.seed,
        check_plots_dir=Path(args.check_plots_dir) if args.check_plots_dir else None,
        no_check_plots=args.no_check_plots,
        extra_header=extra_header,
    )
