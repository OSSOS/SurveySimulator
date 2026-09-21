"""Pure helpers for JWST Sample A grid-cell debiasing (ac2c72)."""
from __future__ import annotations

import math

import numpy as np

A_STEP = 0.2
Q_STEP = 0.2
SI_STEP = 0.001
H_STEP = 0.1
H_COLOR_OFFSET = 1.0
BOWELL_G = -0.12
MOSAIC_AREA_DEG2 = 0.05
MOSAIC_SIDE_DEG = math.sqrt(MOSAIC_AREA_DEG2)
FILL_FACTOR = 1.0
OBLIQUITY_J2000_DEG = 23.4392911
# rot.f95 equat_ecl; used when matching Detos1 / RADECeclXV
F95_OBLIQUITY_ARCSEC = 84381.41
# Eduardo et al. 2026 ICRS mosaic centre (13:57:33, −10:51:55).
# CADC proposal-1568 detector centroids average ~3″ east of this.
FIELD_RA_DEG = 209.3875
FIELD_DEC_DEG = -10.865278
# Orbit-fit reference in Eduardo et al. 2026 §V (JD TDB). This is the
# midpoint of the 10-day campaign, not an observation time.
PAPER_REFERENCE_JD = 2459974.5
# CADC TAP (JWST collection, proposal 1568, NIRCam F150W2 science):
# three 20-tile mosaics, each ~20 h of dithered visits that were
# shift-and-stacked. Use the visit-window midpoint, not 00:00 integer JD.
#   epoch 1  jw01568001*  2023-01-24 09:51 – 01-25 05:34 UTC
#   epoch 2  jw01568002*  2023-01-28 23:38 – 01-29 22:50 UTC
#   epoch 3  jw01568003*  2023-02-04 00:01 – 02-04 19:36 UTC
EPOCH_JD = (2459969.32118, 2459973.96785, 2459979.90854)
# Implant speed range used for characterization (Eduardo et al. 2026 §III.2).
RATE_CUT_MIN_ARCSEC_HR = 0.03
RATE_CUT_MAX_ARCSEC_HR = 8.66
TWO_HOURS_DAY = 2.0 / 24.0


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


def apparent_to_Hr(m_f150w2: float, d_au: float, robs_au: float = 1.0,
                   g: float = BOWELL_G) -> float:
    """H_r inverted from AppMag with r = Δ = d_bary.

    m_r = m_F150W2 + 1 (efficiency-file mapping). Geometry uses the same Bowell
    G=-0.12 law as OSSSSim rather than a constant +0.35 mag phase offset.
    """
    m_r = m_f150w2 + H_COLOR_OFFSET
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


def aimed_at_field(ra_deg: float = FIELD_RA_DEG, dec_deg: float = FIELD_DEC_DEG
                   ) -> tuple[float, float, float, float]:
    """(i, Ω, ω, M) that places a circular orbit on the given ICRS pointing.

    This is the barycentric sky direction, not the apparent direction from JWST.
    Use los_circular_elements to plant in the mosaic as Detos1 sees it.
    """
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
                   obs_jd: float, field_ra: float = FIELD_RA_DEG,
                   field_dec: float = FIELD_DEC_DEG
                   ) -> tuple[float, float, float, float]:
    """Apparent (RA, Dec, sep_deg, rate_arcsec_hr) at obs_jd.

    Detos1 advances M from the element epoch to the pointing JD, then
    measures rate over the next two hours (GetSurvey's second ObsPos).
    """
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
