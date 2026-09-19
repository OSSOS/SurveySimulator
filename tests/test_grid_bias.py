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
        # β ≈ i sin(ω+M); aim puts the object at max |latitude|.
        self.assertAlmostEqual(math.sin(math.radians(arglat)), 1.0 if lat >= 0 else -1.0, places=6)
        lam = (node + arglat) % 360.0
        self.assertAlmostEqual(lam, lon % 360.0, places=4)


if __name__ == "__main__":
    unittest.main()
