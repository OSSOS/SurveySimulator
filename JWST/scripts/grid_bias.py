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
FIELD_RA_DEG = 209.3875
FIELD_DEC_DEG = -10.865278


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

    Uses the small-i approximation λ ≈ Ω + ω + M and β = i sin(ω+M).
    """
    lon, lat = icrs_to_ecliptic(ra_deg, dec_deg)
    inc = max(abs(lat), 0.05)
    arglat = 90.0 if lat >= 0.0 else 270.0
    node = (lon - arglat) % 360.0
    return inc, node, arglat, 0.0


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
