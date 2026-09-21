"""Unit tests for JWST Sample A grid-bias helpers. No Fortran required."""
from __future__ import annotations

import importlib.util
import math
import sys
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
HELPER = ROOT / "JWST" / "scripts" / "grid_bias.py"

spec = importlib.util.spec_from_file_location("grid_bias", HELPER)
grid_bias = importlib.util.module_from_spec(spec)
sys.modules["grid_bias"] = grid_bias
spec.loader.exec_module(grid_bias)


class GridBiasHelpers(unittest.TestCase):
    def test_ifree_zero_on_laplace_plane(self):
        a = 44.0
        ip = grid_bias.laplace_inclination(a)
        om = grid_bias.laplace_node(a)
        self.assertAlmostEqual(grid_bias.compute_ifree(ip, om, a), 0.0, places=6)

    def test_ecliptic_from_ifree_roundtrip(self):
        rng = np.random.default_rng(0)
        a = 44.0
        for ifree in (0.5, 2.0, 5.0, 15.0, 30.0):
            for _ in range(20):
                i, node = grid_bias.ecliptic_from_ifree(ifree, a, rng)
                got = grid_bias.compute_ifree(i, node, a)
                self.assertAlmostEqual(got, ifree, places=5)

    def test_omega_zero_is_not_aligned_with_laplace_node(self):
        a = 43.7
        i = 1.85
        ifree_omega0 = grid_bias.compute_ifree(i, 0.0, a)
        ifree_aligned = grid_bias.compute_ifree(i, grid_bias.laplace_node(a), a)
        self.assertGreater(ifree_omega0, 2.0)
        self.assertLess(ifree_aligned, 0.5)

    def test_hr_inverts_appmag_without_constant_offset(self):
        m_f150w2 = 26.0
        d = 44.0
        h = grid_bias.apparent_to_Hr(m_f150w2, d)
        m_r = m_f150w2 + grid_bias.H_COLOR_OFFSET
        opposition_approx = m_r - 10.0 * math.log10(d)
        # Bowell Φ < 1 at the implied phase, so H is brighter than 10log10(d).
        self.assertLess(h, opposition_approx)
        self.assertGreater(opposition_approx - h, 0.1)
        self.assertLess(opposition_approx - h, 0.5)
        # The old +0.35 term moved H in the opposite direction.
        old = opposition_approx + 0.35
        self.assertGreater(old - h, 0.4)

    def test_sample_aq_near_circular_cell(self):
        rng = np.random.default_rng(1)
        # JPB13-like: a=42.9, e=0, cell a in [42.8, 43.0), q in [42.8, 43.0)
        a, q = grid_bias.sample_aq(rng, (42.8, 43.0), (42.8, 43.0))
        self.assertGreater(a, q)
        self.assertGreater(q, 0.0)
        self.assertGreaterEqual(a, 42.8)
        self.assertLess(a, 43.0)

    def test_sample_aq_rejects_q_greater_than_a_everywhere(self):
        rng = np.random.default_rng(2)
        with self.assertRaises(RuntimeError):
            grid_bias.sample_aq(rng, (40.0, 40.2), (41.0, 41.2), max_tries=50)

    def test_old_098_a_rule_empties_circular_cells(self):
        a0, q0, q1 = 42.8, 42.8, 43.0
        a = a0
        high = min(q1, a * 0.98)
        self.assertLess(high, q0)

    def test_mosaic_fill_is_not_implant_ratio(self):
        self.assertAlmostEqual(grid_bias.MOSAIC_SIDE_DEG ** 2, 0.05, places=12)
        self.assertEqual(grid_bias.FILL_FACTOR, 1.0)
        implant_ff = 0.05 / (1.6 * 1.6)
        self.assertGreater(grid_bias.FILL_FACTOR / implant_ff, 50)

    def test_cell_key_stable(self):
        key = grid_bias.cell_key(44.25, 42.55, 0.0452, 11.37)
        self.assertEqual(key[0], 44.2)
        self.assertEqual(key[1], 42.4)
        bounds = grid_bias.bounds_from_key(key)
        self.assertAlmostEqual(bounds["a"][0], 44.2, places=6)
        self.assertAlmostEqual(bounds["a"][1], 44.4, places=6)
        self.assertAlmostEqual(bounds["Hx"][1] - bounds["Hx"][0], grid_bias.H_STEP)

    def test_figure20_is_not_the_1e5_detection_rate(self):
        lon, lat = grid_bias.icrs_to_ecliptic(
            grid_bias.FIELD_RA_DEG, grid_bias.FIELD_DEC_DEG
        )
        self.assertAlmostEqual(lon, 211.15, places=1)
        self.assertAlmostEqual(lat, 1.07, places=2)
        p7 = grid_bias.geometric_detection_prob(
            grid_bias.MOSAIC_AREA_DEG2, 7.0, lat
        )
        # ~1 per 10^5 draws is the on-sky geometry of 0.05 deg², not Fig. 20.
        self.assertGreater(p7, 5e-6)
        self.assertLess(p7, 2e-5)
        p_cold = grid_bias.geometric_detection_prob(
            grid_bias.MOSAIC_AREA_DEG2, 2.5, lat
        )
        self.assertGreater(p_cold, p7)
        self.assertEqual(
            grid_bias.geometric_detection_prob(
                grid_bias.MOSAIC_AREA_DEG2, 0.5, lat
            ),
            0.0,
        )

    def test_aimed_at_field_reaches_jwst_latitude(self):
        inc, node, peri, M = grid_bias.aimed_at_field()
        lon, lat = grid_bias.icrs_to_ecliptic(
            grid_bias.FIELD_RA_DEG, grid_bias.FIELD_DEC_DEG
        )
        self.assertAlmostEqual(inc, abs(lat), places=5)
        arglat = peri + M
        self.assertAlmostEqual(math.sin(math.radians(arglat)), 1.0 if lat >= 0 else -1.0, places=6)
        lam = (node + arglat) % 360.0
        self.assertAlmostEqual(lam, lon % 360.0, places=4)

    def test_ecliptic_icrf_roundtrip(self):
        x, y, z = -0.55, 0.83, -0.002
        xe, ye, ze = grid_bias.icrf_to_ecliptic(*grid_bias.ecliptic_to_icrf(x, y, z))
        self.assertAlmostEqual(xe, x, places=12)
        self.assertAlmostEqual(ye, y, places=12)
        self.assertAlmostEqual(ze, z, places=12)
        _xi, _yi, zi = grid_bias.ecliptic_to_icrf(x, y, z)
        # JWST near the ecliptic; ICRF z is ~sin(ε) * y ≈ 0.33 AU.
        self.assertGreater(zi, 0.3)

    def test_jwst_csv_observer_is_ecliptic_and_converts(self):
        path = ROOT / "JWST" / "characterization" / "epoch1" / "JWST.csv"
        if not path.is_file():
            self.skipTest(f"missing {path}")
        obs = grid_bias.parse_jpl_horizons_icrf(path, grid_bias.EPOCH_JD[0])
        self.assertAlmostEqual(math.sqrt(sum(c * c for c in obs)), 1.0, places=2)
        self.assertGreater(obs[2], 0.3)

    def test_los_plant_sits_on_field_line_of_sight(self):
        path = ROOT / "JWST" / "characterization" / "epoch1" / "JWST.csv"
        if not path.is_file():
            self.skipTest(f"missing {path}")
        a, e, inc, node, peri, M = grid_bias.los_circular_elements(
            grid_bias.FIELD_RA_DEG, grid_bias.FIELD_DEC_DEG, 44.0, path,
            grid_bias.EPOCH_JD[0],
        )
        self.assertEqual(e, 0.0)
        self.assertAlmostEqual(a, 44.0, places=5)
        # Reconstruct ecliptic position from the circular-element recipe.
        lat = math.degrees(math.asin(math.sin(math.radians(inc))
                                     * math.sin(math.radians(peri + M))))
        lon = (node + peri + M) % 360.0
        x = a * math.cos(math.radians(lon)) * math.cos(math.radians(lat))
        y = a * math.sin(math.radians(lon)) * math.cos(math.radians(lat))
        z = a * math.sin(math.radians(lat))
        obj_icrf = grid_bias.ecliptic_to_icrf(x, y, z)
        obs = grid_bias.parse_jpl_horizons_icrf(path, grid_bias.EPOCH_JD[0])
        los = (
            obj_icrf[0] - obs[0],
            obj_icrf[1] - obs[1],
            obj_icrf[2] - obs[2],
        )
        nrm = math.sqrt(sum(c * c for c in los))
        los = tuple(c / nrm for c in los)
        ra = math.radians(grid_bias.FIELD_RA_DEG)
        dec = math.radians(grid_bias.FIELD_DEC_DEG)
        want = (
            math.cos(dec) * math.cos(ra),
            math.cos(dec) * math.sin(ra),
            math.sin(dec),
        )
        dot = sum(u * v for u, v in zip(los, want))
        sep_deg = math.degrees(math.acos(max(-1.0, min(1.0, dot))))
        self.assertLess(sep_deg, 0.02)

    def test_plant_apparent_radec_is_mosaic_centre(self):
        path = ROOT / "JWST" / "characterization" / "epoch1" / "JWST.csv"
        if not path.is_file():
            self.skipTest(f"missing {path}")
        jd = grid_bias.EPOCH_JD[0]
        a, e, inc, node, peri, M = grid_bias.los_circular_elements(
            grid_bias.FIELD_RA_DEG, grid_bias.FIELD_DEC_DEG, 44.0, path, jd
        )
        obs = grid_bias.parse_jpl_horizons_icrf(path, jd)
        ra, dec = grid_bias.apparent_radec_deg(a, e, inc, node, peri, M, obs)
        sep = grid_bias.sky_separation_deg(
            ra, dec, grid_bias.FIELD_RA_DEG, grid_bias.FIELD_DEC_DEG
        )
        self.assertLess(sep * 60.0, 0.1)
        # Mixed frames (object ICRS, observatory left ecliptic) miss by ~26',
        # larger than the mosaic half-width ~6.7'.
        obs_ecl = grid_bias.icrf_to_ecliptic(*obs)
        ra_m, dec_m = grid_bias.apparent_radec_deg(
            a, e, inc, node, peri, M, obs_ecl
        )
        sep_m = grid_bias.sky_separation_deg(
            ra_m, dec_m, grid_bias.FIELD_RA_DEG, grid_bias.FIELD_DEC_DEG
        )
        self.assertGreater(sep_m, 0.2)
        self.assertGreater(sep_m, grid_bias.MOSAIC_SIDE_DEG / 2.0)

    def test_rate_cut_209_was_field_ra_not_motion_pa(self):
        # Debug log: epoch1 object PA −168.7°, .eff centre 209.4°, hwidth 180°.
        obj = -168.68213509828240
        self.assertGreater(abs(209.4 - obj), 180.0)
        self.assertTrue(grid_bias.angle_in_rate_cone(obj, 209.4, 180.0))
        self.assertTrue(grid_bias.angle_in_rate_cone(obj, 0.0, 180.0))
        # After turnaround the PA is near the old centre, so the unwrapped
        # test passed epochs 2–3 by accident.
        self.assertLess(abs(209.4 - 191.3), 180.0)

    def test_paper_reference_jd_is_not_an_observation(self):
        # Eduardo et al. 2026 §V: JD 2459974.5 is the orbit-fit origin,
        # "near the midpoint of the observation period", not epoch 2.
        self.assertEqual(grid_bias.PAPER_REFERENCE_JD, 2459974.5)
        self.assertNotIn(grid_bias.PAPER_REFERENCE_JD, grid_bias.EPOCH_JD)
        # CADC visit midpoints: ~4.6 d then ~5.9 d, not 1 d and not 5+4 at 00:00.
        d12 = grid_bias.EPOCH_JD[1] - grid_bias.EPOCH_JD[0]
        d23 = grid_bias.EPOCH_JD[2] - grid_bias.EPOCH_JD[1]
        self.assertGreater(d12, 4.0)
        self.assertLess(d12, 5.5)
        self.assertGreater(d23, 5.0)
        self.assertLess(d23, 7.0)
        self.assertAlmostEqual(grid_bias.EPOCH_JD[0], 2459969.32118, places=4)
        self.assertAlmostEqual(grid_bias.EPOCH_JD[2], 2459979.90854, places=4)

    def test_keplerian_plant_stays_in_mosaic_at_cadc_epochs(self):
        path = ROOT / "JWST" / "characterization" / "epoch1" / "JWST.csv"
        if not path.is_file():
            self.skipTest(f"missing {path}")
        element_jd = grid_bias.EPOCH_JD[0]
        a, e, inc, node, peri, M = grid_bias.los_circular_elements(
            grid_bias.FIELD_RA_DEG, grid_bias.FIELD_DEC_DEG, 44.0, path, element_jd
        )
        half = grid_bias.MOSAIC_SIDE_DEG / 2.0
        for jd in grid_bias.EPOCH_JD:
            ra, dec, sep, rate = grid_bias.epoch_geometry(
                a, e, inc, node, peri, M, path, element_jd, jd
            )
            self.assertLess(sep, half, msg=f"JD {jd} sep={sep * 60:.2f}'")
            self.assertGreaterEqual(rate, grid_bias.RATE_CUT_MIN_ARCSEC_HR)
            self.assertLessEqual(rate, grid_bias.RATE_CUT_MAX_ARCSEC_HR)


if __name__ == "__main__":
    unittest.main()
