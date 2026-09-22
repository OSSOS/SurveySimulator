#!/usr/bin/env python3
"""RA / Dec / heliocentric distance of SSim model TNOs near the Roman GBTDS pointing.

This does **not** run Detos1.  Keplerian elements are advanced from the model
epoch to the requested date with the same two-body mean motion Detos1 uses
(n = 360° / (a^{3/2} yr)), then sky coordinates are computed with the Fortran
``pos_cart`` + ``RADECeclXV`` geometry (ecliptic object, ICRF observer).

Default pointing is the Roman GBTDS field given for the JWST follow-up
proposal: 13:52:25.52 −11:01:25.3, epoch 2027-05-01, 10° radius.
"""
from __future__ import annotations

import argparse
import math
from pathlib import Path

import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord, get_body_barycentric
from astropy.table import Table
from astropy.time import Time

# rot.f95 equat_ecl J2000 obliquity
F95_OBLIQUITY_ARCSEC = 84381.41
TWO_PI = 2.0 * math.pi
DEG2RAD = math.pi / 180.0

# Roman GBTDS pointing for the JWST follow-up proposal
GBTDS_RA_HMS = "13:52:25.52"
GBTDS_DEC_DMS = "-11:01:25.3"
DEFAULT_EPOCH = "2027-05-01"
DEFAULT_RADIUS_DEG = 10.0
DEFAULT_MODEL = Path(__file__).resolve().parents[2] / "F95" / "tests" / "Models" / "L7model-3.0-9.0"
REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUT = Path(__file__).resolve().parents[1] / "data" / "gbtds_radec_helio_2027may01.ecsv"


def gbtds_center() -> SkyCoord:
    return SkyCoord(GBTDS_RA_HMS, GBTDS_DEC_DMS, unit=(u.hourangle, u.deg), frame="icrs")


def mean_motion_deg_per_day(a_au: np.ndarray) -> np.ndarray:
    """n = 360° / P, P = a^{3/2} yr in days. Matches Detos1 with gmb ≈ 1."""
    return 360.0 / (np.asarray(a_au, dtype=float) ** 1.5 * 365.25)


def compute_E(e: np.ndarray, M: np.ndarray) -> np.ndarray:
    """Eccentric anomaly via accelerated Newton (Danby / F95 pos_cart)."""
    e = np.asarray(e, dtype=float)
    M = np.mod(np.asarray(M, dtype=float), TWO_PI)
    E = M + 0.85 * np.sign(np.sin(M)) * e
    for _ in range(20):
        sin_e = e * np.sin(E)
        f = E - sin_e - M
        iterate = np.abs(f) > 1e-14
        if not np.any(iterate):
            break
        cos_e = e[iterate] * np.cos(E[iterate])
        fp = 1.0 - cos_e
        fpp = sin_e[iterate]
        fppp = cos_e
        de = -f[iterate] / fp
        de = -f[iterate] / (fp + de * fpp / 2.0)
        de = -f[iterate] / (fp + de * fpp / 2.0 + de * de * fppp / 6.0)
        E = E.copy()
        E[iterate] = E[iterate] + de
    else:
        raise ValueError("POS_CART: eccentric anomaly did not converge")
    return E


def pos_cart(a: np.ndarray, e: np.ndarray, inc: np.ndarray,
             node: np.ndarray, peri: np.ndarray, M: np.ndarray
             ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Heliocentric ecliptic J2000 (x, y, z) AU. Signs match F95/elemutils.f95."""
    a = np.asarray(a, dtype=float)
    e = np.asarray(e, dtype=float)
    inc = np.asarray(inc, dtype=float)
    node = np.asarray(node, dtype=float)
    peri = np.asarray(peri, dtype=float)
    M = np.mod(np.asarray(M, dtype=float), TWO_PI)
    signe = np.sign(a)
    signe = np.where(signe == 0.0, 1.0, signe)
    cos_i = np.cos(inc)
    sin_i = np.sqrt(np.clip(1.0 - cos_i * cos_i, 0.0, None))
    c_w, s_w = np.cos(peri), np.sin(peri)
    c_o, s_o = np.cos(node), np.sin(node)
    delau6 = signe * np.sqrt(a * signe)
    delau7 = np.abs(delau6) * np.sqrt((1.0 - e * e) * signe)
    E = compute_E(e, M)
    q0 = delau6 ** 2 * (np.cos(E) - e)
    q1 = delau7 * delau6 * np.sin(E)
    x = (c_o * c_w - cos_i * s_o * s_w) * q0 + (-c_o * s_w - cos_i * s_o * c_w) * q1
    y = (s_o * c_w + cos_i * c_o * s_w) * q0 + (-s_o * s_w + cos_i * c_o * c_w) * q1
    z = (sin_i * s_w) * q0 + (sin_i * c_w) * q1
    return x, y, z


def ecliptic_to_icrf(x: np.ndarray, y: np.ndarray, z: np.ndarray
                     ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """equat_ecl(-1): conventional ecliptic J2000 → ICRF."""
    eps = math.radians(F95_OBLIQUITY_ARCSEC / 3600.0)
    coseps, sineps = math.cos(eps), math.sin(eps)
    return x, coseps * y - sineps * z, sineps * y + coseps * z


def helio_dist_au(x: np.ndarray, y: np.ndarray, z: np.ndarray) -> np.ndarray:
    return np.sqrt(x * x + y * y + z * z)


def apparent_radec_deg(x_ecl: np.ndarray, y_ecl: np.ndarray, z_ecl: np.ndarray,
                       obs_icrf: np.ndarray
                       ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """ICRS RA/Dec and geocentric distance as RADECeclXV computes them."""
    ox, oy, oz = ecliptic_to_icrf(x_ecl, y_ecl, z_ecl)
    relx = ox - obs_icrf[0]
    rely = oy - obs_icrf[1]
    relz = oz - obs_icrf[2]
    delta = np.sqrt(relx * relx + rely * rely + relz * relz)
    ra = np.degrees(np.arctan2(rely, relx)) % 360.0
    dec = np.degrees(np.arcsin(np.clip(relz / delta, -1.0, 1.0)))
    return ra, dec, delta


def sky_separation_deg(ra1: np.ndarray, dec1: np.ndarray,
                       ra2: float, dec2: float) -> np.ndarray:
    c1 = SkyCoord(ra=ra1 * u.deg, dec=dec1 * u.deg, frame="icrs")
    c2 = SkyCoord(ra=ra2 * u.deg, dec=dec2 * u.deg, frame="icrs")
    return c1.separation(c2).to(u.deg).value


def observer_icrf_au(epoch: Time) -> np.ndarray:
    """Geocenter, heliocentric ICRF AU (Sun at origin, matching Keplerian pos_cart)."""
    earth = get_body_barycentric("earth", epoch)
    sun = get_body_barycentric("sun", epoch)
    return (earth - sun).xyz.to(u.au).value


def parse_model_epoch_jd(path: Path) -> float:
    with path.open() as f_obj:
        for line in f_obj:
            if not line.startswith("#"):
                break
            if "JD" in line and "=" in line:
                return float(line.split("=")[1].split()[0].replace("d", "e"))
    raise ValueError(f"No JD epoch in header of {path}")


def load_ssim_model(path: Path) -> dict:
    """Load an L7-style SSim model (comment header + whitespace columns)."""
    a, e, inc, node, peri, M, H, dist, comp, j, k = (
        [], [], [], [], [], [], [], [], [], [], []
    )
    with path.open() as f_obj:
        for line in f_obj:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split()
            if len(parts) < 8:
                continue
            a.append(float(parts[0]))
            e.append(float(parts[1]))
            inc.append(float(parts[2]))
            node.append(float(parts[3]))
            peri.append(float(parts[4]))
            M.append(float(parts[5]))
            H.append(float(parts[6]))
            dist.append(float(parts[7]))
            comp.append(parts[8] if len(parts) > 8 else "unknown")
            j.append(int(parts[9]) if len(parts) > 9 else 0)
            k.append(int(parts[10]) if len(parts) > 10 else 0)
    return {
        "a": np.asarray(a),
        "e": np.asarray(e),
        "inc": np.asarray(inc),
        "node": np.asarray(node),
        "peri": np.asarray(peri),
        "M": np.asarray(M),
        "H": np.asarray(H),
        "dist_model": np.asarray(dist),
        "comp": np.asarray(comp),
        "j": np.asarray(j, dtype=int),
        "k": np.asarray(k, dtype=int),
        "epoch_jd": parse_model_epoch_jd(path),
        "filename": str(path),
        "n_model": len(a),
    }


def positions_at_epoch(model: dict, obs_jd: float, obs_icrf: np.ndarray) -> dict:
    """Advance M and return apparent RA/Dec plus heliocentric r at obs_jd."""
    dt_day = obs_jd - model["epoch_jd"]
    m_obs = model["M"] + mean_motion_deg_per_day(model["a"]) * dt_day
    x, y, z = pos_cart(
        model["a"],
        model["e"],
        model["inc"] * DEG2RAD,
        model["node"] * DEG2RAD,
        model["peri"] * DEG2RAD,
        m_obs * DEG2RAD,
    )
    r = helio_dist_au(x, y, z)
    ra, dec, delta = apparent_radec_deg(x, y, z, obs_icrf)
    out = dict(model)
    out.update({
        "M_obs": np.mod(m_obs, 360.0),
        "x": x,
        "y": y,
        "z": z,
        "helio_dist": r,
        "ra": ra,
        "dec": dec,
        "delta": delta,
        "obs_jd": obs_jd,
    })
    return out


def select_field(positions: dict, ra0: float, dec0: float, radius_deg: float
                 ) -> dict:
    sep = sky_separation_deg(positions["ra"], positions["dec"], ra0, dec0)
    keep = sep <= radius_deg
    out = {
        key: (val[keep] if isinstance(val, np.ndarray) and val.shape[:1] == (len(sep),)
              else val)
        for key, val in positions.items()
    }
    out["sep_deg"] = sep[keep]
    out["n_in_field"] = int(keep.sum())
    out["n_model"] = positions["n_model"]
    out["radius_deg"] = radius_deg
    out["field_ra"] = ra0
    out["field_dec"] = dec0
    return out


def _repo_relative(path: str) -> str:
    try:
        return str(Path(path).resolve().relative_to(REPO_ROOT))
    except ValueError:
        return path


def field_table(selected: dict) -> Table:
    table = Table()
    table["ra"] = selected["ra"]
    table["dec"] = selected["dec"]
    table["helio_dist"] = selected["helio_dist"]
    table["delta"] = selected["delta"]
    table["sep_deg"] = selected["sep_deg"]
    table["a"] = selected["a"]
    table["e"] = selected["e"]
    table["inc"] = selected["inc"]
    table["node"] = selected["node"]
    table["peri"] = selected["peri"]
    table["M"] = selected["M"]
    table["M_obs"] = selected["M_obs"]
    table["H"] = selected["H"]
    table["comp"] = selected["comp"]
    table["j"] = selected["j"]
    table["k"] = selected["k"]
    table["ra"].unit = u.deg
    table["dec"].unit = u.deg
    table["helio_dist"].unit = u.au
    table["delta"].unit = u.au
    table["sep_deg"].unit = u.deg
    table["a"].unit = u.au
    table["inc"].unit = u.deg
    table["node"].unit = u.deg
    table["peri"].unit = u.deg
    table["M"].unit = u.deg
    table["M_obs"].unit = u.deg
    table["H"].unit = u.mag
    table.meta = {
        "description": "SSim L7 RA/Dec/heliocentric distance near the Roman GBTDS pointing",
        "field_ra_hms": GBTDS_RA_HMS,
        "field_dec_dms": GBTDS_DEC_DMS,
        "field_ra_deg": selected["field_ra"],
        "field_dec_deg": selected["field_dec"],
        "radius_deg": selected["radius_deg"],
        "obs_jd": selected["obs_jd"],
        "obs_iso": Time(selected["obs_jd"], format="jd").isot,
        "model_epoch_jd": selected["epoch_jd"],
        "model_file": _repo_relative(selected["filename"]),
        "n_model": selected["n_model"],
        "n_in_field": selected["n_in_field"],
        "observer": "geocenter heliocentric ICRF",
        "geometry": "SSim pos_cart + RADECeclXV; no Detos1 / detectability",
    }
    return table


def write_plots(table: Table, out_prefix: Path, ra0: float, dec0: float,
                radius_deg: float) -> list[Path]:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    ra = np.asarray(table["ra"], dtype=float)
    dec = np.asarray(table["dec"], dtype=float)
    r = np.asarray(table["helio_dist"], dtype=float)
    paths = []

    fig, ax = plt.subplots(figsize=(8.0, 7.2))
    sc = ax.scatter(ra, dec, c=r, s=8, cmap="viridis", alpha=0.85, linewidths=0)
    ax.add_patch(plt.Circle((ra0, dec0), radius_deg, fill=False, color="k",
                            lw=1.2, ls="--"))
    ax.plot(ra0, dec0, marker="+", color="k", ms=10, mew=1.4)
    ax.set_xlabel("RA [deg, ICRS]")
    ax.set_ylabel("Dec [deg, ICRS]")
    ax.set_title("SSim L7 near Roman GBTDS  2027-05-01")
    ax.set_aspect("equal", adjustable="box")
    ax.invert_xaxis()
    cb = fig.colorbar(sc, ax=ax, pad=0.02)
    cb.set_label("heliocentric distance [AU]")
    ax.legend(handles=[Line2D([0], [0], color="k", ls="--", lw=1.2,
                              label=f"{radius_deg:.0f}° radius")],
              loc="upper right", fontsize=9)
    fig.tight_layout()
    sky_path = Path(str(out_prefix) + "_sky.png")
    fig.savefig(sky_path, dpi=140)
    plt.close(fig)
    paths.append(sky_path)

    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    ax.hist(r, bins=min(40, max(12, int(np.sqrt(len(r))))), color="0.35",
            histtype="stepfilled", alpha=0.85)
    ax.set_xlabel("heliocentric distance [AU]")
    ax.set_ylabel("N")
    ax.set_title(f"helio_dist in {radius_deg:.0f}° of GBTDS  (n={len(r)})")
    fig.tight_layout()
    dist_path = Path(str(out_prefix) + "_helio_dist.png")
    fig.savefig(dist_path, dpi=140)
    plt.close(fig)
    paths.append(dist_path)
    return paths


def summarize(table: Table) -> str:
    ra = np.asarray(table["ra"], dtype=float)
    dec = np.asarray(table["dec"], dtype=float)
    r = np.asarray(table["helio_dist"], dtype=float)
    comps, counts = np.unique(table["comp"], return_counts=True)
    lines = [
        f"n_in_field = {len(table)}  /  n_model = {table.meta['n_model']}",
        f"epoch      = {table.meta['obs_iso']}  (JD {table.meta['obs_jd']:.5f})",
        f"pointing   = {table.meta['field_ra_hms']} {table.meta['field_dec_dms']}"
        f"  ({table.meta['field_ra_deg']:.5f}, {table.meta['field_dec_deg']:.5f})",
        f"radius     = {table.meta['radius_deg']} deg",
        f"RA         = {ra.min():.4f} … {ra.max():.4f} deg",
        f"Dec        = {dec.min():.4f} … {dec.max():.4f} deg",
        f"helio_dist = {r.min():.2f} … {r.max():.2f} AU  (median {np.median(r):.2f})",
        "components:",
    ]
    for name, n in zip(comps, counts):
        lines.append(f"  {name:12s}  {n:6d}")
    return "\n".join(lines)


def parse_args(argv=None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--model", type=Path, default=DEFAULT_MODEL,
                   help="SSim model file (L7-style)")
    p.add_argument("--epoch", default=DEFAULT_EPOCH,
                   help="observation epoch (ISO date or JD)")
    p.add_argument("--radius", type=float, default=DEFAULT_RADIUS_DEG,
                   help="field radius in degrees")
    p.add_argument("--ra", default=GBTDS_RA_HMS,
                   help="field RA (hms or degrees)")
    p.add_argument("--dec", default=GBTDS_DEC_DMS,
                   help="field Dec (dms or degrees)")
    p.add_argument("--out", type=Path, default=DEFAULT_OUT,
                   help="output ECSV path")
    p.add_argument("--no-plots", action="store_true")
    return p.parse_args(argv)


def _parse_skycoord(ra_str: str, dec_str: str) -> SkyCoord:
    try:
        return SkyCoord(ra_str, dec_str, unit=(u.hourangle, u.deg), frame="icrs")
    except Exception:
        return SkyCoord(float(ra_str) * u.deg, float(dec_str) * u.deg, frame="icrs")


def _parse_epoch(value: str) -> Time:
    try:
        jd = float(value)
        if jd > 2000000:
            return Time(jd, format="jd")
    except ValueError:
        pass
    return Time(value)


def main(argv=None) -> int:
    args = parse_args(argv)
    field = _parse_skycoord(args.ra, args.dec)
    epoch = _parse_epoch(args.epoch)
    model = load_ssim_model(args.model)
    obs_icrf = observer_icrf_au(epoch)
    positions = positions_at_epoch(model, epoch.jd, obs_icrf)
    selected = select_field(
        positions, field.ra.deg, field.dec.deg, args.radius
    )
    table = field_table(selected)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    table.write(args.out, format="ascii.ecsv", overwrite=True)
    csv_path = args.out.with_suffix(".csv")
    table["ra", "dec", "helio_dist"].write(csv_path, format="ascii.csv", overwrite=True)
    print(summarize(table))
    print(f"wrote {args.out}")
    print(f"wrote {csv_path}")
    if not args.no_plots:
        prefix = args.out.with_suffix("")
        for path in write_plots(table, prefix, field.ra.deg, field.dec.deg, args.radius):
            print(f"wrote {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
