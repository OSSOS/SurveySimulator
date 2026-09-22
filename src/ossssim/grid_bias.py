"""Grid-cell Horvitz–Thompson helpers shared by JWST Sample A and N26.

Module-level field defaults are JWST Sample A so existing JWST scripts and
tests keep working.  N26 (and any later survey) must pass ``survey=`` or
explicit field geometry.
"""
from __future__ import annotations

import csv
import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np

A_STEP = 0.2
Q_STEP = 0.2
SI_STEP = 0.001
H_STEP = 0.1
BOWELL_G = -0.12
OBLIQUITY_J2000_DEG = 23.4392911
# rot.f95 equat_ecl; used when matching Detos1 / RADECeclXV
F95_OBLIQUITY_ARCSEC = 84381.41
POINTINGS_TEMPLATE_NAME = "pointings.template"
TWO_HOURS_DAY = 2.0 / 24.0
# Napier et al. 2026 / user: r_AB = STMAG_F606W - 0.3
STMAG_F606W_TO_R_AB = -0.3


@dataclass(frozen=True)
class GridSurvey:
    """Pencil-beam survey geometry and photometry mapping for grid debiasing."""

    name: str
    field_ra_deg: float
    field_dec_deg: float
    mosaic_width_deg: float
    mosaic_height_deg: float
    epoch_jd: tuple
    mag_color_offset: float
    mag_column: str
    observer_csv: str
    eff_file: str
    fill_factor: float = 1.0
    paper_reference_jd: float | None = None
    rate_cut_min_arcsec_hr: float = 0.03
    rate_cut_max_arcsec_hr: float = 8.66
    epoch_layout: str = "subdir"  # "subdir" → epoch{i}/; "flat" → char root
    detections_relpath: str = "data/detections.csv"
    detections_full_name: str = "detections-full"
    check_detected_title: str = "detected flag≥4"

    @property
    def n_epochs(self) -> int:
        return len(self.epoch_jd)

    @property
    def mosaic_area_deg2(self) -> float:
        return self.mosaic_width_deg * self.mosaic_height_deg

    @property
    def mosaic_side_deg(self) -> float:
        """Side of the equal-area square; JWST mosaic is actually square."""
        return math.sqrt(self.mosaic_area_deg2)

    def mag_to_r(self, mag: float) -> float:
        """Map the survey catalog magnitude onto the OSSOS r_AB system."""
        return mag + self.mag_color_offset


# Eduardo et al. 2026 ICRS mosaic centre (13:57:33, −10:51:55).
# CADC proposal-1568 detector centroids average ~3″ east of this.
JWST_SAMPLE_A = GridSurvey(
    name="JWST Sample A",
    field_ra_deg=209.3875,
    field_dec_deg=-10.865278,
    mosaic_width_deg=math.sqrt(0.05),
    mosaic_height_deg=math.sqrt(0.05),
    epoch_jd=(2459969.32118, 2459973.96785, 2459979.90854),
    mag_color_offset=1.0,  # m_r = m_F150W2 + 1
    mag_column="m_f150w2",
    observer_csv="JWST.csv",
    eff_file="JWST_sampleA.eff",
    paper_reference_jd=2459974.5,
    rate_cut_min_arcsec_hr=0.03,
    rate_cut_max_arcsec_hr=8.66,
    epoch_layout="subdir",
    detections_relpath="data/jwst_sampleA.csv",
    detections_full_name="JWST-free-cla_m.detections-full",
    check_detected_title="detected flag≥4 at all 3 epochs",
)

# Napier et al. 2026 (PSJ 7, 117) reanalysis of Bernstein et al. 2004
# GO-9433 ACS WFC F606W. 6-tile 400″×600″ mosaic, ~0.02 deg².
# Single 15-day heliostack at the full-span midpoint, not 3-epoch AND.
N26_HELIOSTACK = GridSurvey(
    name="N26 heliostack",
    field_ra_deg=15.0 * (14.0 + 7.0 / 60.0 + 53.33 / 3600.0),
    field_dec_deg=-(11.0 + 21.0 / 60.0 + 38.0 / 3600.0),
    mosaic_width_deg=400.0 / 3600.0,
    mosaic_height_deg=600.0 / 3600.0,
    epoch_jd=(2452672.8585,),
    mag_color_offset=STMAG_F606W_TO_R_AB,  # r_AB = STMAG_F606W - 0.3
    mag_column="m_stmag",
    observer_csv="HST.csv",
    eff_file="N26.eff",
    paper_reference_jd=2452672.8585,
    rate_cut_min_arcsec_hr=0.05,
    rate_cut_max_arcsec_hr=6.4,
    epoch_layout="flat",
    detections_relpath="data/n26_detections.csv",
    detections_full_name="N26-free-cla_m.detections-full",
    check_detected_title="detected flag≥4 (single 15-day stack)",
)

# JWST aliases kept so JWST/scripts/grid_bias.py can re-export this module.
H_COLOR_OFFSET = JWST_SAMPLE_A.mag_color_offset
MOSAIC_AREA_DEG2 = JWST_SAMPLE_A.mosaic_area_deg2
MOSAIC_SIDE_DEG = JWST_SAMPLE_A.mosaic_side_deg
MOSAIC_WIDTH_DEG = JWST_SAMPLE_A.mosaic_width_deg
MOSAIC_HEIGHT_DEG = JWST_SAMPLE_A.mosaic_height_deg
FILL_FACTOR = JWST_SAMPLE_A.fill_factor
FIELD_RA_DEG = JWST_SAMPLE_A.field_ra_deg
FIELD_DEC_DEG = JWST_SAMPLE_A.field_dec_deg
PAPER_REFERENCE_JD = JWST_SAMPLE_A.paper_reference_jd
EPOCH_JD = JWST_SAMPLE_A.epoch_jd
RATE_CUT_MIN_ARCSEC_HR = JWST_SAMPLE_A.rate_cut_min_arcsec_hr
RATE_CUT_MAX_ARCSEC_HR = JWST_SAMPLE_A.rate_cut_max_arcsec_hr


def stmag_f606w_to_r_ab(stmag: float) -> float:
    """AB r from ACS F606W STMAG. Napier ST plus the −0.3 zeropoint."""
    return float(stmag) + STMAG_F606W_TO_R_AB


def _survey_or_default(survey: GridSurvey | None) -> GridSurvey:
    return JWST_SAMPLE_A if survey is None else survey


def laplace_inclination(a_au: float) -> float:
    return 1.759 + 0.0321 * (a_au - 41.8)


def laplace_node(a_au: float) -> float:
    return 90.0 - 0.5 * (a_au - 43.0)


def compute_ifree(i_deg: float, node_deg: float, a_au: float) -> float:
    """Free inclination relative to the Laplace plane, given ecliptic (i, Ω)."""
    ip = laplace_inclination(a_au)
    om_lp = laplace_node(a_au)
    cos_ifree = (
        math.cos(math.radians(i_deg)) * math.cos(math.radians(ip))
        + math.sin(math.radians(i_deg)) * math.sin(math.radians(ip))
        * math.cos(math.radians(node_deg - om_lp))
    )
    return math.degrees(math.acos(max(-1.0, min(1.0, cos_ifree))))


def _orbit_pole(i_deg: float, node_deg: float) -> np.ndarray:
    i = math.radians(i_deg)
    node = math.radians(node_deg)
    return np.array([
        math.sin(i) * math.sin(node),
        -math.sin(i) * math.cos(node),
        math.cos(i),
    ])


def ecliptic_from_ifree(ifree_deg: float, a_au: float, rng: np.random.Generator) -> tuple[float, float]:
    """Draw ecliptic (i, Ω) whose free inclination is ifree_deg.

    The orbit pole is placed at angular distance i_free from the Laplace pole
    with a uniform random azimuth. Inverse of compute_ifree.
    """
    ip = laplace_inclination(a_au)
    om_lp = laplace_node(a_au)
    lp = _orbit_pole(ip, om_lp)
    ref = np.array([1.0, 0.0, 0.0]) if abs(lp[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
    e1 = np.cross(lp, ref)
    e1 /= np.linalg.norm(e1)
    e2 = np.cross(lp, e1)
    psi = rng.uniform(0.0, 2.0 * math.pi)
    ifree = math.radians(ifree_deg)
    pole = math.cos(ifree) * lp + math.sin(ifree) * (math.cos(psi) * e1 + math.sin(psi) * e2)
    pole = pole / np.linalg.norm(pole)
    i_ecl = math.degrees(math.acos(max(-1.0, min(1.0, float(pole[2])))))
    si = math.sin(math.radians(i_ecl))
    if si < 1e-12:
        node_ecl = float(rng.uniform(0.0, 360.0))
    else:
        node_ecl = math.degrees(math.atan2(pole[0] / si, -pole[1] / si)) % 360.0
    return max(i_ecl, 0.05), node_ecl


def bowell_phase_correction(alpha_rad: float, g: float = BOWELL_G) -> float:
    """2.5 log10((1-G)φ1 + G φ2), the AppMag term added when inverting for H."""
    if alpha_rad <= 0.0:
        return 0.0
    ta = math.tan(alpha_rad / 2.0)
    phi1 = math.exp(-3.33 * ta ** 0.63)
    phi2 = math.exp(-1.87 * ta ** 1.22)
    phi = (1.0 - g) * phi1 + g * phi2
    if phi <= 0.0:
        return 0.0
    return 2.5 * math.log10(phi)


def apparent_to_Hr(m_survey: float, d_au: float, robs_au: float = 1.0,
                   g: float = BOWELL_G, color_offset: float | None = None,
                   survey: GridSurvey | None = None) -> float:
    """H_r inverted from AppMag with r = Δ = d_bary.

    ``color_offset`` maps the survey magnitude onto r_AB (JWST: +1 for
    F150W2; N26: −0.3 for STMAG F606W). Geometry uses the same Bowell
    G=-0.12 law as OSSSSim rather than a constant +0.35 mag phase offset.
    """
    if color_offset is None:
        color_offset = _survey_or_default(survey).mag_color_offset
    m_r = m_survey + color_offset
    denom = 2.0 * d_au * d_au
    cos_a = max(-1.0, min(1.0, (-robs_au ** 2 + 2.0 * d_au ** 2) / denom))
    alpha = math.acos(cos_a)
    return m_r - 5.0 * math.log10(d_au * d_au) + bowell_phase_correction(alpha, g)


def geometric_detection_prob(area_deg2: float, inc_deg: float, beta_deg: float) -> float:
    """Single-epoch geometric probability for a small field.

    P ≈ A / (360° × 2 × sqrt(i² − β²)) when i > |β|. For the JWST mosaic
    (A=0.05 deg²) and a ~7° inclination belt this is ~1e-5. Eduardo et al.
    2026 Figure 20 is the cold-belt H_r luminosity function, not this rate.
    """
    if abs(inc_deg) <= abs(beta_deg):
        return 0.0
    return area_deg2 / (360.0 * 2.0 * math.sqrt(inc_deg ** 2 - beta_deg ** 2))


def geometric_prob_for_aimed(a: float, e: float, inc_deg: float, node_deg: float,
                             peri_deg: float, M_deg: float,
                             area_deg2: float | None = None,
                             survey: GridSurvey | None = None) -> float:
    """P_geom at the aimed orbit's own ecliptic latitude, not the ICRS field.

    Parallax (~1°) means barycentric β of an on-LOS TNO is not icrs_to_ecliptic
    of the pointing. Using the field latitude spuriously zeros P_geom for cold
    cells where aimed i sits between |β_obj| and |β_field|.
    """
    if area_deg2 is None:
        area_deg2 = _survey_or_default(survey).mosaic_area_deg2
    x, y, z = ecliptic_xyz_from_elements(a, e, inc_deg, node_deg, peri_deg, M_deg)
    r = math.sqrt(x * x + y * y + z * z)
    if r <= 0.0:
        return 0.0
    beta = math.degrees(math.asin(max(-1.0, min(1.0, z / r))))
    return geometric_detection_prob(area_deg2, inc_deg, beta)


def icrs_to_ecliptic(ra_deg: float, dec_deg: float) -> tuple[float, float]:
    """J2000 equatorial (RA, Dec) to ecliptic (lon, lat), degrees."""
    ra = math.radians(ra_deg)
    dec = math.radians(dec_deg)
    eps = math.radians(OBLIQUITY_J2000_DEG)
    x = math.cos(dec) * math.cos(ra)
    y = math.cos(dec) * math.sin(ra)
    z = math.sin(dec)
    ye = y * math.cos(eps) + z * math.sin(eps)
    ze = -y * math.sin(eps) + z * math.cos(eps)
    lon = math.degrees(math.atan2(ye, x)) % 360.0
    lat = math.degrees(math.asin(max(-1.0, min(1.0, ze))))
    return lon, lat


def aimed_at_field(ra_deg: float | None = None, dec_deg: float | None = None,
                   survey: GridSurvey | None = None
                   ) -> tuple[float, float, float, float]:
    """(i, Ω, ω, M) that places a circular orbit on the given ICRS pointing.

    This is the barycentric sky direction, not the apparent direction from JWST.
    Use los_circular_elements to plant in the mosaic as Detos1 sees it.
    """
    surv = _survey_or_default(survey)
    if ra_deg is None:
        ra_deg = surv.field_ra_deg
    if dec_deg is None:
        dec_deg = surv.field_dec_deg
    lon, lat = icrs_to_ecliptic(ra_deg, dec_deg)
    inc = max(abs(lat), 0.05)
    arglat = 90.0 if lat >= 0.0 else 270.0
    node = (lon - arglat) % 360.0
    return inc, node, arglat, 0.0


def _obliquity_rad() -> float:
    return math.radians(F95_OBLIQUITY_ARCSEC / 3600.0)


def ecliptic_to_icrf(x: float, y: float, z: float) -> tuple[float, float, float]:
    """equat_ecl(-1): ecliptic J2000 → ICRF."""
    coseps = math.cos(_obliquity_rad())
    sineps = math.sin(_obliquity_rad())
    return x, coseps * y - sineps * z, sineps * y + coseps * z


def icrf_to_ecliptic(x: float, y: float, z: float) -> tuple[float, float, float]:
    """equat_ecl(+1): ICRF → ecliptic J2000."""
    coseps = math.cos(_obliquity_rad())
    sineps = math.sin(_obliquity_rad())
    return x, coseps * y + sineps * z, -sineps * y + coseps * z


def circular_elements_through_ecliptic_xyz(
        x: float, y: float, z: float) -> tuple[float, float, float, float, float, float]:
    """Circular (a, e, i, Ω, ω, M) whose position is ecliptic (x,y,z) AU."""
    r = math.sqrt(x * x + y * y + z * z)
    lat = math.degrees(math.asin(max(-1.0, min(1.0, z / r))))
    lon = math.degrees(math.atan2(y, x)) % 360.0
    inc = max(abs(lat), 0.05)
    arglat = 90.0 if lat >= 0.0 else 270.0
    node = (lon - arglat) % 360.0
    return r, 0.0, inc, node, arglat, 0.0


def parse_jpl_horizons_icrf(path, jd: float) -> tuple[float, float, float]:
    """Observer barycentric ICRF (AU) from a Horizons CSV, matching read_jpl_csv."""
    with open(path) as fh:
        text = fh.read()
    header, _, rest = text.partition("$$SOE")
    ecliptic_frame = False
    for line in header.splitlines():
        if "Reference frame" in line:
            ecliptic_frame = ("Ecliptic" in line) or ("ecliptic" in line)
    for line in rest.splitlines():
        raw = line.strip()
        if not raw:
            continue
        if raw.startswith("$$EOE"):
            break
        parts = [p.strip() for p in raw.split(",")]
        if len(parts) < 8:
            continue
        ejd = float(parts[0])
        if ejd <= jd:
            continue
        x, y, z = (float(parts[i]) for i in (2, 3, 4))
        vx, vy, vz = (float(parts[i]) for i in (5, 6, 7))
        dt = jd - ejd
        pos = (x + vx * dt, y + vy * dt, z + vz * dt)
        if ecliptic_frame:
            pos = ecliptic_to_icrf(*pos)
        return pos
    raise RuntimeError(f"JD {jd} outside Horizons range in {path}")


def los_circular_elements(ra_deg: float, dec_deg: float, a_au: float,
                          jpl_path, jd: float
                          ) -> tuple[float, float, float, float, float, float]:
    """Circular orbit at a_au along the observer LOS to (ra, dec) at jd."""
    obs = parse_jpl_horizons_icrf(jpl_path, jd)
    ra = math.radians(ra_deg)
    dec = math.radians(dec_deg)
    los = (
        math.cos(dec) * math.cos(ra),
        math.cos(dec) * math.sin(ra),
        math.sin(dec),
    )
    b = 2.0 * sum(o * l for o, l in zip(obs, los))
    c = sum(o * o for o in obs) - a_au * a_au
    disc = max(0.0, b * b - 4.0 * c)
    t = 0.5 * (-b + math.sqrt(disc))
    obj_icrf = tuple(o + t * l for o, l in zip(obs, los))
    obj_ecl = icrf_to_ecliptic(*obj_icrf)
    return circular_elements_through_ecliptic_xyz(*obj_ecl)


def ecliptic_xyz_from_elements(a: float, e: float, inc_deg: float,
                               node_deg: float, peri_deg: float, M_deg: float
                               ) -> tuple[float, float, float]:
    """Barycentric ecliptic xyz matching F95 pos_cart (e near 0 is fine)."""
    inc = math.radians(inc_deg)
    node = math.radians(node_deg)
    peri = math.radians(peri_deg)
    M = math.radians(M_deg) % (2.0 * math.pi)
    E = M
    for _ in range(20):
        f = E - e * math.sin(E) - M
        if abs(f) < 1e-14:
            break
        E -= f / (1.0 - e * math.cos(E))
    cos_i, sin_i = math.cos(inc), math.sin(inc)
    c_w, s_w = math.cos(peri), math.sin(peri)
    c_o, s_o = math.cos(node), math.sin(node)
    q0 = a * (math.cos(E) - e)
    q1 = a * math.sqrt(max(0.0, 1.0 - e * e)) * math.sin(E)
    x = (c_o * c_w - cos_i * s_o * s_w) * q0 + (-c_o * s_w - cos_i * s_o * c_w) * q1
    y = (s_o * c_w + cos_i * c_o * s_w) * q0 + (-s_o * s_w + cos_i * c_o * c_w) * q1
    z = (sin_i * s_w) * q0 + (sin_i * c_w) * q1
    return x, y, z


def apparent_radec_deg(a: float, e: float, inc_deg: float, node_deg: float,
                       peri_deg: float, M_deg: float, obs_icrf
                       ) -> tuple[float, float]:
    """ICRS RA/Dec as RADECeclXV computes them (object ecliptic, obs ICRF)."""
    obj_icrf = ecliptic_to_icrf(*ecliptic_xyz_from_elements(
        a, e, inc_deg, node_deg, peri_deg, M_deg))
    rel = [obj_icrf[i] - obs_icrf[i] for i in range(3)]
    delta = math.sqrt(sum(v * v for v in rel))
    ra = math.degrees(math.atan2(rel[1], rel[0])) % 360.0
    dec = math.degrees(math.asin(max(-1.0, min(1.0, rel[2] / delta))))
    return ra, dec


def sky_separation_deg(ra1: float, dec1: float, ra2: float, dec2: float) -> float:
    dra = (ra1 - ra2) * math.cos(math.radians(0.5 * (dec1 + dec2)))
    return math.hypot(dra, dec1 - dec2)


def angle_in_rate_cone(obj_deg: float, centre_deg: float, hwidth_deg: float
                       ) -> bool:
    """Whether a motion PA is inside the rate_cut direction cone.

    Detos1 uses atan2 ∈ [−180°, 180°]. A centre of 209.4° (the JWST field
    RA, not a PA) compared without wrapping rejects pre-turnaround motion
    at −168.7° even when half-width is 180°.
    """
    dang = (centre_deg - obj_deg + 180.0) % 360.0 - 180.0
    return abs(dang) <= hwidth_deg


def mean_motion_deg_per_day(a_au: float) -> float:
    """n = 360° / P, P = a^{3/2} yr in days. Matches Detos1 with gmb≈1."""
    return 360.0 / (a_au ** 1.5 * 365.25)


def epoch_geometry(a: float, e: float, inc_deg: float, node_deg: float,
                   peri_deg: float, M_deg: float, jpl_path, element_jd: float,
                   obs_jd: float, field_ra: float | None = None,
                   field_dec: float | None = None,
                   survey: GridSurvey | None = None
                   ) -> tuple[float, float, float, float]:
    """Apparent (RA, Dec, sep_deg, rate_arcsec_hr) at obs_jd.

    Detos1 advances M from the element epoch to the pointing JD, then
    measures rate over the next two hours (GetSurvey's second ObsPos).
    """
    surv = _survey_or_default(survey)
    if field_ra is None:
        field_ra = surv.field_ra_deg
    if field_dec is None:
        field_dec = surv.field_dec_deg
    n = mean_motion_deg_per_day(a)
    m1 = M_deg + n * (obs_jd - element_jd)
    m2 = M_deg + n * (obs_jd + TWO_HOURS_DAY - element_jd)
    obs1 = parse_jpl_horizons_icrf(jpl_path, obs_jd)
    obs2 = parse_jpl_horizons_icrf(jpl_path, obs_jd + TWO_HOURS_DAY)
    ra1, dec1 = apparent_radec_deg(a, e, inc_deg, node_deg, peri_deg, m1, obs1)
    ra2, dec2 = apparent_radec_deg(a, e, inc_deg, node_deg, peri_deg, m2, obs2)
    dra = (ra1 - ra2) * math.cos(math.radians(dec1))
    ddec = dec2 - dec1
    rate = math.hypot(dra, ddec) / TWO_HOURS_DAY * 3600.0 / 24.0
    sep = sky_separation_deg(ra1, dec1, field_ra, field_dec)
    return ra1, dec1, sep, rate


def cell_index(value: float, step: float) -> float:
    return math.floor(value / step) * step


def cell_key(a: float, q: float, sin_ifree: float, hx: float) -> tuple:
    return (
        round(cell_index(a, A_STEP), 6),
        round(cell_index(q, Q_STEP), 6),
        round(cell_index(sin_ifree, SI_STEP), 6),
        round(cell_index(hx, H_STEP), 6),
    )


def bounds_from_key(key: tuple) -> dict:
    a0, q0, si0, h0 = key
    return {
        "a": (a0, a0 + A_STEP),
        "q": (q0, q0 + Q_STEP),
        "sin_ifree": (max(0.0, si0), min(1.0, si0 + SI_STEP)),
        "Hx": (h0, h0 + H_STEP),
    }


def load_detections(path, survey: GridSurvey) -> list[dict]:
    """Read a survey detections CSV and assign (a, q, sin i_free, H) cells.

    Uses ``survey.mag_column`` and ``survey.mag_color_offset`` for H_r.
    If the CSV has an ``ifree`` column (N26 Table 2 midpoints), that value
    is used; otherwise i_free is computed from ecliptic i with Ω=0.
    """
    rows = []
    with Path(path).open() as fh:
        for row in csv.DictReader(fh):
            a, e, i = float(row["a"]), float(row["e"]), float(row["i"])
            d = float(row["d_bary"])
            mag = float(row[survey.mag_column])
            hx = apparent_to_Hr(mag, d, survey=survey)
            q = a * (1.0 - e)
            if "ifree" in row and str(row["ifree"]).strip():
                ifree = float(row["ifree"])
            else:
                ifree = compute_ifree(i, 0.0, a)
            rows.append({
                **row, "a": a, "e": e, "i": i, "d_bary": d, "q": q,
                "Hx": hx, "ifree": ifree, "mag": mag,
                "sin_ifree": math.sin(math.radians(ifree)),
                "cell": cell_key(a, q, math.sin(math.radians(ifree)), hx),
            })
    return rows


def sample_aq(rng: np.random.Generator, a_bounds: tuple, q_bounds: tuple,
              max_tries: int = 10000) -> tuple[float, float]:
    """Uniform draw in the (a, q) rectangle restricted to 0 < q < a."""
    a0, a1 = a_bounds
    q0, q1 = q_bounds
    for _ in range(max_tries):
        a = float(rng.uniform(a0, a1))
        q = float(rng.uniform(q0, q1))
        if 0.0 < q < a:
            return a, q
    raise RuntimeError(
        f"empty (a,q) cell a=[{a0}, {a1}) q=[{q0}, {q1}); no bound orbit with q < a"
    )


def render_pointings_text(template: str, epoch: int, jd: float,
                          survey: GridSurvey | None = None) -> str:
    """Fill pointings.template for one epoch. GetSurvey reads pointings.list."""
    surv = _survey_or_default(survey)
    text = template.format(
        epoch=epoch,
        jd=jd,
        ra=surv.field_ra_deg,
        dec=surv.field_dec_deg,
        side=surv.mosaic_side_deg,
        width=surv.mosaic_width_deg,
        height=surv.mosaic_height_deg,
        fill=surv.fill_factor,
        observer_csv=surv.observer_csv,
        eff_file=surv.eff_file,
    )
    if not text.endswith("\n"):
        text += "\n"
    return text


def setup_pointings(char_root, template_path=None,
                    survey: GridSurvey | None = None) -> list:
    """Write gitignored pointings.list from characterization/pointings.template.

    Detos1/GetSurvey always open `{survey_dir}/pointings.list`. Keep the
    committed source as a template and regenerate the list at run time so
    pulling this branch does not require resetting those files.

    ``epoch_layout='subdir'`` (JWST) writes ``epoch{i}/pointings.list``.
    ``epoch_layout='flat'`` (N26) writes ``pointings.list`` in ``char_root``.
    """
    surv = _survey_or_default(survey)
    char_root = Path(char_root)
    template_path = Path(template_path) if template_path else (
        char_root / POINTINGS_TEMPLATE_NAME
    )
    template = template_path.read_text()
    written = []
    if surv.epoch_layout == "flat":
        dests = [(char_root / "pointings.list", 1, surv.epoch_jd[0])]
    else:
        dests = [
            (char_root / f"epoch{idx}" / "pointings.list", idx, jd)
            for idx, jd in enumerate(surv.epoch_jd, start=1)
        ]
    for dest, idx, jd in dests:
        dest.parent.mkdir(parents=True, exist_ok=True)
        text = render_pointings_text(template, idx, jd, survey=surv)
        if not dest.exists() or dest.read_text() != text:
            dest.write_text(text)
        written.append(dest)
    return written


def icrs_los_unit(ra_deg: float, dec_deg: float) -> np.ndarray:
    """ICRS unit vector at (RA, Dec). Detos1 FoV tests are in this frame."""
    ra = math.radians(ra_deg)
    dec = math.radians(dec_deg)
    return np.array([
        math.cos(dec) * math.cos(ra),
        math.cos(dec) * math.sin(ra),
        math.sin(dec),
    ])


def barycentric_on_icrs_los(obs_icrf, ra_deg: float, dec_deg: float, r_au: float
                            ) -> np.ndarray | None:
    """Far |R|=r intersection of the ICRS LOS, returned in J2000 ecliptic.

    Observatory vectors from JWST.csv are converted to ICRF; the FoV (RA, Dec)
    is ICRS. Orbit elements are ecliptic, so the intersection is rotated with
    icrf_to_ecliptic (the same trick as los_circular_elements).
    """
    los = icrs_los_unit(ra_deg, dec_deg)
    obs = np.asarray(obs_icrf, dtype=float)
    b = 2.0 * float(obs @ los)
    c = float(obs @ obs) - r_au * r_au
    disc = b * b - 4.0 * c
    if disc < 0.0:
        return None
    root = math.sqrt(disc)
    t = 0.5 * (-b + root)
    if t <= 0.0:
        t = 0.5 * (-b - root)
        if t <= 0.0:
            return None
    obj_icrf = obs + t * los
    return np.array(icrf_to_ecliptic(*obj_icrf))


def sample_orbital_radius(a: float, e: float, rng: np.random.Generator) -> float:
    """Draw barycentric r on the ellipse, r ∈ [q, Q] = [a(1-e), a(1+e)]."""
    if e < 1e-12:
        return a
    q = a * (1.0 - e)
    q_ap = a * (1.0 + e)
    if q_ap <= q:
        return a
    return float(rng.uniform(q, q_ap))


def true_anomaly_from_radius(a: float, e: float, r_au: float) -> float:
    """|f| in radians from the orbit equation. Caller chooses the sign of f."""
    if e < 1e-12:
        return 0.0
    cos_f = (a * (1.0 - e * e) / r_au - 1.0) / e
    return math.acos(max(-1.0, min(1.0, cos_f)))


def mean_anomaly_from_true(e: float, f_rad: float) -> float:
    """M = E − e sin E, with E from true anomaly f (radians)."""
    if e < 1e-12:
        return f_rad
    cos_f = math.cos(f_rad)
    sin_f = math.sin(f_rad)
    den = 1.0 + e * cos_f
    cos_E = (e + cos_f) / den
    sin_E = math.sqrt(max(0.0, 1.0 - e * e)) * sin_f / den
    ecc = math.atan2(sin_E, cos_E)
    return ecc - e * math.sin(ecc)


def argument_of_latitude(R, inc_deg: float, node_deg: float) -> float:
    """u = ω + f in radians from ecliptic position and (i, Ω)."""
    inc = math.radians(inc_deg)
    node = math.radians(node_deg)
    xhat = np.array([math.cos(node), math.sin(node), 0.0])
    yhat = np.array([
        -math.cos(inc) * math.sin(node),
        math.cos(inc) * math.cos(node),
        math.sin(inc),
    ])
    R = np.asarray(R, dtype=float)
    return math.atan2(float(R @ yhat), float(R @ xhat))


def _radius_on_ellipse(a: float, e: float, r_au: float) -> float | None:
    q = a * (1.0 - e)
    q_ap = a * (1.0 + e)
    if r_au < q - 1e-8 or r_au > q_ap + 1e-8:
        return None
    return min(q_ap, max(q, r_au))


def nodes_from_inclination(R, inc_deg: float) -> list[float]:
    """Ascending nodes Ω (deg) whose plane of inclination i contains R.

    n · R = 0 with n = (sin i sin Ω, −sin i cos Ω, cos i). Empty if |β| > i.
    """
    x, y, z = (float(c) for c in np.asarray(R, dtype=float))
    r = math.hypot(math.hypot(x, y), z)
    if r < 1e-18:
        return []
    inc = math.radians(inc_deg)
    si, ci = math.sin(inc), math.cos(inc)
    if abs(si) < 1e-12:
        return [0.0] if abs(z / r) < 1e-8 else []
    # x sin Ω − y cos Ω = −z cot i
    amp_a, amp_b, target = x, -y, -z * ci / si
    amp = math.hypot(amp_a, amp_b)
    if amp < 1e-18 or abs(target) > amp + 1e-10:
        return []
    psi = math.atan2(amp_b, amp_a)
    alpha = math.asin(max(-1.0, min(1.0, target / amp)))
    out = []
    seen = set()
    for ang in (alpha - psi, math.pi - alpha - psi):
        node = math.degrees(ang) % 360.0
        key = round(node, 8)
        if key in seen:
            continue
        seen.add(key)
        out.append(node)
    return out


def peri_m_from_position(a: float, e: float, inc_deg: float, node_deg: float,
                         R, f_sign: float = 1.0) -> tuple[float, float]:
    """(ω, M) in degrees from (a, e, i, Ω) and ecliptic position R.

    r = |R| fixes |true anomaly| f; u is the argument of latitude of R;
    ω = u − f; M follows from f.
    """
    r_au = float(np.linalg.norm(R))
    f_abs = true_anomaly_from_radius(a, e, r_au)
    f_rad = f_abs if f_sign >= 0.0 else -f_abs
    u_rad = argument_of_latitude(R, inc_deg, node_deg)
    peri = math.degrees(u_rad - f_rad) % 360.0
    mean_anom = math.degrees(mean_anomaly_from_true(e, f_rad)) % 360.0
    return peri, mean_anom


def keplerian_at_radec_r(a: float, e: float, inc_deg: float,
                         ra_deg: float, dec_deg: float, r_au: float,
                         obs_icrf, f_sign: float = 1.0, node_index: int = 0
                         ) -> tuple[float, float, float, float, float, float] | None:
    """Map (a, e, i) and ICRS (RA, Dec, r) to (a, e, i, Ω, ω, M).

    This is the reusable geometric step. H is not used (photometry only).
    RA/Dec are ICRS; elements are J2000 ecliptic. The LOS is intersected
    at |R|=r in ICRF and rotated with icrf_to_ecliptic.

    Returns None if r is off the ellipse, the ray misses the sphere, or
    |β| > i so no node exists. Two Ω solutions in general; node_index
    selects one. f_sign chooses inbound vs outbound true anomaly.
    """
    r_au = _radius_on_ellipse(a, e, r_au)
    if r_au is None:
        return None
    pos_ecl = barycentric_on_icrs_los(obs_icrf, ra_deg, dec_deg, r_au)
    if pos_ecl is None:
        return None
    nodes = nodes_from_inclination(pos_ecl, inc_deg)
    if not nodes:
        return None
    node = nodes[int(node_index) % len(nodes)]
    peri, mean_anom = peri_m_from_position(
        a, e, inc_deg, node, pos_ecl, f_sign
    )
    return a, e, inc_deg, node, peri, mean_anom


def poles_through_position_at_ifree(R, ifree_deg: float, a_au: float) -> list:
    """Orbit poles P with P·R = 0 and angle(P, Laplace pole) = i_free.

    Two solutions in general (the plane can tilt either way around R). Empty
    if |β_lp| > i_free, so the line of sight cannot sit on that free
    inclination.
    """
    lp = _orbit_pole(laplace_inclination(a_au), laplace_node(a_au))
    rhat = np.asarray(R, dtype=float)
    nrm = np.linalg.norm(rhat)
    if nrm < 1e-18:
        return []
    rhat = rhat / nrm
    ref = np.array([1.0, 0.0, 0.0]) if abs(rhat[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
    e1 = np.cross(rhat, ref)
    n1 = np.linalg.norm(e1)
    if n1 < 1e-18:
        return []
    e1 = e1 / n1
    e2 = np.cross(rhat, e1)
    amp_a = float(e1 @ lp)
    amp_b = float(e2 @ lp)
    amp = math.hypot(amp_a, amp_b)
    target = math.cos(math.radians(ifree_deg))
    if amp < 1e-18 or abs(target) > amp + 1e-10:
        return []
    phi = math.atan2(amp_b, amp_a)
    dth = math.acos(max(-1.0, min(1.0, target / amp)))
    out = []
    seen = []
    for th in (phi + dth, phi - dth):
        pole = math.cos(th) * e1 + math.sin(th) * e2
        pole = pole / np.linalg.norm(pole)
        key = tuple(round(float(c), 10) for c in pole)
        if key in seen:
            continue
        seen.append(key)
        i_ecl = math.degrees(math.acos(max(-1.0, min(1.0, float(pole[2])))))
        si = math.sin(math.radians(i_ecl))
        if si < 1e-12:
            node_ecl = 0.0
        else:
            node_ecl = math.degrees(math.atan2(pole[0] / si, -pole[1] / si)) % 360.0
        out.append((max(i_ecl, 0.05), node_ecl))
    return out


def sample_mosaic_icrs(rng: np.random.Generator,
                       ra_deg: float | None = None,
                       dec_deg: float | None = None,
                       side_deg: float | None = None,
                       width_deg: float | None = None,
                       height_deg: float | None = None,
                       survey: GridSurvey | None = None
                       ) -> tuple[float, float]:
    """Uniform ICRS (RA, Dec) inside the rectangular mosaic pointing."""
    surv = _survey_or_default(survey)
    if ra_deg is None:
        ra_deg = surv.field_ra_deg
    if dec_deg is None:
        dec_deg = surv.field_dec_deg
    if width_deg is None:
        width_deg = side_deg if side_deg is not None else surv.mosaic_width_deg
    if height_deg is None:
        height_deg = side_deg if side_deg is not None else surv.mosaic_height_deg
    return (
        float(ra_deg + rng.uniform(-0.5 * width_deg, 0.5 * width_deg)),
        float(dec_deg + rng.uniform(-0.5 * height_deg, 0.5 * height_deg)),
    )


def aimed_elements(a: float, e: float, ifree_deg: float,
                   ra_deg: float, dec_deg: float, r_au: float,
                   obs_icrf, f_sign: float = 1.0, pole_index: int = 0
                   ) -> tuple[float, float, float, float] | None:
    """Solve (i, Ω, ω, M) for a grid-cell i_free at ICRS (RA, Dec, r).

    i_free plus the ecliptic position fixes (i, Ω) (the orbit pole through
    R at i_free from the Laplace pole). Then peri_m_from_position gives
    (ω, M). Use keplerian_at_radec_r when ecliptic i is already known.
    """
    r_au = _radius_on_ellipse(a, e, r_au)
    if r_au is None:
        return None
    pos_ecl = barycentric_on_icrs_los(obs_icrf, ra_deg, dec_deg, r_au)
    if pos_ecl is None:
        return None
    poles = poles_through_position_at_ifree(pos_ecl, ifree_deg, a)
    if not poles:
        return None
    inc, node = poles[int(pole_index) % len(poles)]
    peri, mean_anom = peri_m_from_position(
        a, e, inc, node, pos_ecl, f_sign
    )
    return inc, node, peri, mean_anom


def sample_aimed_elements(a: float, e: float, ifree_deg: float, obs_icrf,
                          rng: np.random.Generator, max_tries: int = 40,
                          survey: GridSurvey | None = None
                          ) -> tuple[float, float, float, float] | None:
    """FoV-aimed (i, Ω, ω, M) for one (a, e, i_free) draw.

    Samples an ICRS location in the mosaic and r on [q, Q], then inverts.
    Discrete branches (±f, two poles) are chosen uniformly. Returns None
    if no inversion succeeds (cheap: i_free below the Laplace latitude).
    """
    for _ in range(max_tries):
        ra, dec = sample_mosaic_icrs(rng, survey=survey)
        r_au = sample_orbital_radius(a, e, rng)
        f_sign = 1.0 if rng.random() < 0.5 else -1.0
        pole_index = int(rng.integers(0, 2))
        el = aimed_elements(
            a, e, ifree_deg, ra, dec, r_au, obs_icrf, f_sign, pole_index
        )
        if el is not None:
            return el
    return None


def aimed_detection_bias(n_aimed: int, geom_weight_sum: float) -> float:
    """Horvitz–Thompson P(detect | cell) from FoV-aimed draws.

    n_detected / n_aimed is P(Sample A | FoV). Each aimed orbit carries
    geometric_detection_prob(A, i, β) so the product is P(detect | cell),
    the same quantity isotropic (Ω, ω, M) sampling estimates ~1e5× slower.
    geom_weight_sum is Σ 1_detected P_geom over aimed draws.
    """
    if n_aimed <= 0:
        return 0.0
    return float(geom_weight_sum) / float(n_aimed)


CHECK_PLOT_KEYS = ("ra", "dec", "a", "e", "i", "Omega", "omega", "M")


def empty_check_samples() -> dict:
    return {key: [] for key in CHECK_PLOT_KEYS}


def record_check_sample(store: dict, ra: float, dec: float, a: float, e: float,
                        inc: float, node: float, peri: float, M: float) -> None:
    store["ra"].append(float(ra))
    store["dec"].append(float(dec))
    store["a"].append(float(a))
    store["e"].append(float(e))
    store["i"].append(float(inc))
    store["Omega"].append(float(node) % 360.0)
    store["omega"].append(float(peri) % 360.0)
    store["M"].append(float(M) % 360.0)


def as_check_arrays(samples: dict) -> dict:
    return {key: np.asarray(samples[key], dtype=float) for key in CHECK_PLOT_KEYS}


def stack_check_samples(parts: list) -> dict:
    """Concatenate per-cell check-sample dicts for a run-level plot."""
    out = {}
    for key in CHECK_PLOT_KEYS:
        chunks = [np.asarray(part[key], dtype=float) for part in parts]
        chunks = [c for c in chunks if c.size]
        out[key] = np.concatenate(chunks) if chunks else np.array([], dtype=float)
    return out


def check_plot_tag(label) -> str:
    """Filename stem for check plots: Sample A object name, or cell key."""
    if isinstance(label, tuple):
        return "cell_" + "_".join(f"{float(v):.4g}" for v in label).replace(".", "p")
    text = str(label).strip() or "object"
    return "".join(c if c.isalnum() or c in "-_" else "_" for c in text)


def write_bias_check_plots(out_dir, sampled: dict, detected: dict, tag: str,
                           field_ra: float | None = None,
                           field_dec: float | None = None,
                           side_deg: float | None = None,
                           width_deg: float | None = None,
                           height_deg: float | None = None,
                           survey: GridSurvey | None = None) -> list:
    """RA/Dec and a/e/i/Ω/ω/M check plots: sampled vs detected (flag≥4).

    Sampled = aimed orbits sent through Detos1. Detected = those with
    flag≥4 (all epochs for a multi-epoch survey). `tag` is the filename
    stem. The two RA/Dec clouds should fill the mosaic the same way if
    detection is not a spatial cut inside the field; the element
    histograms should match if detection is not a function of those
    elements.
    """
    import matplotlib
    matplotlib.use("Agg")
    from matplotlib import pyplot as plt
    from matplotlib.patches import Rectangle
    from pathlib import Path

    surv = _survey_or_default(survey)
    if field_ra is None:
        field_ra = surv.field_ra_deg
    if field_dec is None:
        field_dec = surv.field_dec_deg
    if width_deg is None:
        width_deg = side_deg if side_deg is not None else surv.mosaic_width_deg
    if height_deg is None:
        height_deg = side_deg if side_deg is not None else surv.mosaic_height_deg
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    sampled = as_check_arrays(sampled)
    detected = as_check_arrays(detected)
    n_s = int(sampled["ra"].size)
    n_d = int(detected["ra"].size)
    half_w = 0.5 * width_deg
    half_h = 0.5 * height_deg
    pad_w = 0.4 * width_deg
    pad_h = 0.4 * height_deg
    ra_lim = (field_ra - half_w - pad_w, field_ra + half_w + pad_w)
    dec_lim = (field_dec - half_h - pad_h, field_dec + half_h + pad_h)

    fig, axes = plt.subplots(2, 2, figsize=(9.2, 8.0),
                             gridspec_kw={"height_ratios": [1.35, 1.0]})
    scat_s, scat_d = axes[0]
    hist_ra, hist_dec = axes[1]
    for ax, data, n, color, title in (
            (scat_s, sampled, n_s, "0.35", f"sampled for Detos1 (n={n_s})"),
            (scat_d, detected, n_d, "C0", f"detected flag≥4 (n={n_d})"),
    ):
        if n:
            ax.scatter(data["ra"], data["dec"], s=6, alpha=0.35, c=color,
                       linewidths=0, rasterized=True)
        ax.add_patch(Rectangle(
            (field_ra - half_w, field_dec - half_h), width_deg, height_deg,
            fill=False, edgecolor="k", lw=1.0,
        ))
        ax.plot(field_ra, field_dec, "k+", ms=9, mew=1.2)
        ax.set_xlim(*ra_lim)
        ax.set_ylim(*dec_lim)
        ax.set_aspect("equal", adjustable="box")
        ax.set_xlabel("RA [deg, ICRS]")
        ax.set_ylabel("Dec [deg, ICRS]")
        ax.set_title(title)
        ax.grid(True, alpha=0.25)
    _overlay_hist(hist_ra, sampled["ra"], detected["ra"], "RA [deg, ICRS]", n_s, n_d)
    _overlay_hist(hist_dec, sampled["dec"], detected["dec"], "Dec [deg, ICRS]", n_s, n_d)
    fig.suptitle(f"{tag}: epoch-1 RA/Dec  ({surv.check_detected_title})",
                 fontsize=11)
    fig.tight_layout()
    radec_path = out_dir / f"check_{tag}_radec.png"
    fig.savefig(radec_path, dpi=140)
    plt.close(fig)

    elem_labels = (
        ("a", "a [au]", False),
        ("e", "e", False),
        ("i", "i [deg]", False),
        ("Omega", r"$\Omega$ [deg]", True),
        ("omega", r"$\omega$ [deg]", True),
        ("M", "M [deg]", True),
    )
    fig, axes = plt.subplots(2, 3, figsize=(10.5, 6.4))
    for ax, (key, xlabel, circular) in zip(axes.ravel(), elem_labels):
        _overlay_hist(ax, sampled[key], detected[key], xlabel, n_s, n_d,
                      circular=circular)
    fig.suptitle(
        f"{tag}: elements sent to Detos1 vs detected (flag≥4). "
        "Densities should match if detection is independent of these elements.",
        fontsize=10,
    )
    fig.tight_layout()
    elem_path = out_dir / f"check_{tag}_elements.png"
    fig.savefig(elem_path, dpi=140)
    plt.close(fig)
    return [radec_path, elem_path]


def _overlay_hist(ax, sampled, detected, xlabel: str, n_s: int, n_d: int,
                  circular: bool = False) -> None:
    sampled = np.asarray(sampled, dtype=float)
    detected = np.asarray(detected, dtype=float)
    if circular:
        hist_range = (0.0, 360.0)
        bins = 36
    elif sampled.size:
        lo, hi = float(np.min(sampled)), float(np.max(sampled))
        if hi <= lo:
            hi = lo + 1e-6
        pad = 0.05 * (hi - lo)
        hist_range = (lo - pad, hi + pad)
        bins = min(40, max(12, int(np.sqrt(sampled.size))))
    else:
        hist_range = None
        bins = 20
    if sampled.size:
        ax.hist(sampled, bins=bins, range=hist_range, density=True,
                histtype="stepfilled", alpha=0.35, color="0.45",
                label=f"sampled ({n_s})")
    if detected.size:
        ax.hist(detected, bins=bins, range=hist_range, density=True,
                histtype="step", color="C0", lw=1.6,
                label=f"detected flag≥4 ({n_d})")
    ax.set_xlabel(xlabel)
    ax.set_ylabel("density")
    ax.legend(fontsize=7, frameon=False)
    if circular:
        ax.set_xlim(0.0, 360.0)
