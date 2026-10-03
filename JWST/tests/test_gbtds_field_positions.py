"""Project tests for JWST/scripts/gbtds_field_positions.py (no Fortran)."""
from __future__ import annotations

import math
import sys
import unittest
from pathlib import Path

import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.time import Time

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "JWST" / "scripts"))
import gbtds_field_positions as gfp  # noqa: E402

TEST_MODEL = ROOT / "tests" / "data" / "test_model.dat"


class PosCartTest(unittest.TestCase):
    def test_perihelion_on_x_axis(self):
        x, y, z = gfp.pos_cart(
            np.array([10.0]), np.array([0.2]), np.array([0.0]),
            np.array([0.0]), np.array([0.0]), np.array([0.0]),
        )
        self.assertAlmostEqual(float(x[0]), 8.0, places=10)
        self.assertAlmostEqual(float(y[0]), 0.0, places=10)
        self.assertAlmostEqual(float(z[0]), 0.0, places=10)

    def test_mean_motion_includes_gmb(self):
        n = gfp.mean_motion_deg_per_day(np.array([39.4]))
        period_day = 39.4 ** 1.5 * 365.25
        self.assertAlmostEqual(float(n[0]) * period_day, 360.0 * math.sqrt(gfp.GMB),
                               places=8)


class FieldSelectTest(unittest.TestCase):
    def test_gbtds_center_parses(self):
        c = gfp.gbtds_center()
        self.assertAlmostEqual(c.ra.deg, 208.1063333, places=5)
        self.assertAlmostEqual(c.dec.deg, -11.0236944, places=5)

    def test_numeric_ra_is_degrees(self):
        c = gfp._parse_skycoord("208.1063333", "-11.0236944")
        self.assertAlmostEqual(c.ra.deg, 208.1063333, places=5)
        self.assertAlmostEqual(c.dec.deg, -11.0236944, places=5)

    def test_radius_cut(self):
        center = gfp.gbtds_center()
        n = 8
        dummy = {
            "ra": np.array([center.ra.deg, center.ra.deg + 20.0] + [center.ra.deg] * (n - 2)),
            "dec": np.array([center.dec.deg] * n),
            "n_model": n,
            "a": np.ones(n),
            "filename": "x",
            "epoch_jd": 2453157.5,
        }
        for key in ("e", "inc", "node", "peri", "M", "H", "dist_model",
                    "helio_dist", "delta", "M_obs", "x", "y", "z"):
            dummy[key] = np.zeros(n)
        dummy["comp"] = np.array(["t"] * n)
        dummy["j"] = np.zeros(n, dtype=int)
        dummy["k"] = np.zeros(n, dtype=int)
        dummy["obs_jd"] = 2461516.5
        sel = gfp.select_field(dummy, center.ra.deg, center.dec.deg, 10.0)
        self.assertEqual(sel["n_in_field"], n - 1)
        self.assertTrue(np.all(sel["sep_deg"] <= 10.0))


class TinyModelEphemerisTest(unittest.TestCase):
    def test_test_model_runs_and_has_finite_coords(self):
        model = gfp.load_ssim_model(TEST_MODEL)
        self.assertEqual(model["n_model"], 2)
        epoch = Time("2027-05-01")
        pos = gfp.positions_at_epoch(model, epoch.jd)
        self.assertEqual(pos["ra"].shape, (2,))
        self.assertTrue(np.all(np.isfinite(pos["ra"])))
        self.assertTrue(np.all(np.isfinite(pos["dec"])))
        self.assertTrue(np.all(pos["helio_dist"] > 1.0))

    def test_object_on_los_lands_near_field(self):
        """Circular orbit along the GBTDS ICRS LOS at barycentric r = 42 AU."""
        epoch = Time("2027-05-01")
        earth = gfp.observer_icrf_au(epoch)
        sun = gfp.barycentric_icrf_au("sun", epoch)
        field = gfp.gbtds_center()
        r = 42.0
        ra, dec = field.ra.radian, field.dec.radian
        los = np.array([math.cos(dec) * math.cos(ra),
                        math.cos(dec) * math.sin(ra),
                        math.sin(dec)])
        b = 2.0 * np.dot(earth, los)
        c = np.dot(earth, earth) - r * r
        t = 0.5 * (-b + math.sqrt(b * b - 4.0 * c))
        obj_icrf = earth + t * los
        helio_expected = float(np.linalg.norm(obj_icrf - sun))
        eps = math.radians(gfp.F95_OBLIQUITY_ARCSEC / 3600.0)
        x, y0, z0 = obj_icrf
        y = math.cos(eps) * y0 + math.sin(eps) * z0
        z = -math.sin(eps) * y0 + math.cos(eps) * z0
        lat = math.degrees(math.asin(z / r))
        lon = math.degrees(math.atan2(y, x)) % 360.0
        inc = max(abs(lat), 0.05)
        arglat = 90.0 if lat >= 0.0 else 270.0
        node = (lon - arglat) % 360.0
        model = {
            "a": np.array([r]),
            "e": np.array([0.0]),
            "inc": np.array([inc]),
            "node": np.array([node]),
            "peri": np.array([arglat]),
            "M": np.array([0.0]),
            "H": np.array([8.0]),
            "dist_model": np.array([r]),
            "comp": np.array(["planted"]),
            "j": np.array([0]),
            "k": np.array([0]),
            "epoch_jd": epoch.jd,
            "filename": "planted",
            "n_model": 1,
        }
        pos = gfp.positions_at_epoch(model, epoch.jd, earth, sun)
        sep = SkyCoord(pos["ra"] * u.deg, pos["dec"] * u.deg).separation(field)
        self.assertLess(sep.deg[0], 0.05)
        self.assertAlmostEqual(float(pos["helio_dist"][0]), helio_expected, places=4)


if __name__ == "__main__":
    unittest.main()
