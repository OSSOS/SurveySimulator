"""Unit tests for Napier et al. 2026 N26 heliostack grid-bias characterization."""
from __future__ import annotations

import csv
import importlib.util
import math
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
N26_ROOT = ROOT / "N26"
HELPER = N26_ROOT / "scripts" / "grid_bias.py"
CHAR_PATH = ROOT / "src" / "ossssim" / "characterization.py"


def _load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


grid_bias = _load("n26_grid_bias", HELPER)
characterization = _load("ossssim_characterization_impl", CHAR_PATH)
N26 = grid_bias.N26_HELIOSTACK
SingleTanhParam = characterization.SingleTanhParam
logistic_to_tanh = characterization.logistic_to_tanh


class NapierPhotometry(unittest.TestCase):
    def test_stmag_to_r_ab_offset(self):
        self.assertEqual(grid_bias.STMAG_F606W_TO_R_AB, -0.3)
        self.assertAlmostEqual(grid_bias.stmag_f606w_to_r_ab(29.21), 28.91)
        self.assertAlmostEqual(N26.mag_to_r(29.21), 28.91)
        self.assertAlmostEqual(N26.mag_color_offset, -0.3)

    def test_survey_m50_is_r_ab(self):
        A, m0, w = logistic_to_tanh(1.0, grid_bias.stmag_f606w_to_r_ab(29.21), 0.11)
        self.assertAlmostEqual(A, 1.0)
        self.assertAlmostEqual(m0, 28.91)
        self.assertAlmostEqual(w, 0.22)

    def test_logistic_matches_ossos_tanh(self):
        eta0, m50, sigma = 1.0, 28.91, 0.11
        A, m0, w = logistic_to_tanh(eta0, m50, sigma)
        mag = np.linspace(27.5, 30.0, 50)
        logistic = eta0 / (1.0 + np.exp((mag - m50) / sigma))
        tanh = 0.5 * A * (1.0 - np.tanh((mag - m0) / w))
        np.testing.assert_allclose(tanh, logistic, rtol=0, atol=1e-12)

    def test_n26_eff_matches_napier_mapped_to_r(self):
        text = (N26_ROOT / "characterization" / "N26.eff").read_text()
        line = [ln for ln in text.splitlines() if ln.startswith("single_param")][0]
        param = SingleTanhParam.from_string(line)
        A, m0, w = logistic_to_tanh(1.0, grid_bias.stmag_f606w_to_r_ab(29.21), 0.11)
        self.assertAlmostEqual(param.A, A, places=2)
        self.assertAlmostEqual(param.M_0, m0, places=2)
        self.assertAlmostEqual(param.w, w, places=2)
        self.assertAlmostEqual(param.efficiency(m0), 0.5, places=6)
        self.assertGreater(param.efficiency(m0 + 0.3), 0.05)
        self.assertLess(param.efficiency(m0 + 0.3), 0.25)

    def test_hr_uses_minus_0p3_not_plus_one(self):
        d = 42.19
        h_n26 = grid_bias.apparent_to_Hr(27.99, d)
        h_if_jwst = grid_bias.apparent_to_Hr(27.99, d, color_offset=1.0, survey=None)
        self.assertLess(h_n26, h_if_jwst - 1.0)
        opposition = 27.69 - 10.0 * math.log10(d)
        self.assertLess(h_n26, opposition)
        self.assertGreater(opposition - h_n26, 0.1)


class NapierField(unittest.TestCase):
    def test_bernstein_pointing_and_area(self):
        ra = N26.field_ra_deg
        dec = N26.field_dec_deg
        self.assertAlmostEqual(ra, 15.0 * (14 + 7 / 60 + 53.33 / 3600), places=6)
        self.assertAlmostEqual(dec, -(11 + 21 / 60 + 38 / 3600), places=6)
        self.assertAlmostEqual(N26.mosaic_width_deg, 400.0 / 3600.0)
        self.assertAlmostEqual(N26.mosaic_height_deg, 600.0 / 3600.0)
        area = N26.mosaic_area_deg2
        self.assertGreater(area, 0.018)
        self.assertLess(area, 0.021)
        self.assertEqual(N26.n_epochs, 1)
        self.assertEqual(N26.epoch_layout, "flat")
        self.assertAlmostEqual(N26.epoch_jd[0], 2452672.8585, places=4)

    def test_geometric_p_is_smaller_than_jwst(self):
        _, lat = grid_bias.icrs_to_ecliptic(N26.field_ra_deg, N26.field_dec_deg)
        p = grid_bias.geometric_detection_prob(N26.mosaic_area_deg2, 7.0, lat)
        self.assertGreater(p, 1e-6)
        self.assertLess(p, 1e-5)

    def test_hst_csv_covers_midpoint_and_rate_dt(self):
        path = N26_ROOT / "characterization" / "HST.csv"
        self.assertTrue(path.is_file())
        jd = N26.epoch_jd[0]
        obs = grid_bias.parse_jpl_horizons_icrf(path, jd)
        r = math.sqrt(sum(c * c for c in obs))
        # Early February: Earth is near perihelion (~0.98 au), not 1.00.
        self.assertGreater(r, 0.97)
        self.assertLess(r, 1.03)
        obs2 = grid_bias.parse_jpl_horizons_icrf(path, jd + 2.0 / 24.0)
        r2 = math.sqrt(sum(c * c for c in obs2))
        self.assertGreater(r2, 0.97)
        self.assertLess(r2, 1.03)

    def test_setup_pointings_is_flat_single_epoch(self):
        src = N26_ROOT / "characterization" / "pointings.template"
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "pointings.template").write_text(src.read_text())
            paths = grid_bias.setup_pointings(root)
            self.assertEqual(len(paths), 1)
            self.assertEqual(paths[0].name, "pointings.list")
            self.assertFalse((root / "epoch1").exists())
            text = paths[0].read_text()
            data = [ln for ln in text.splitlines() if ln and not ln.startswith("#")]
            self.assertEqual(len(data), 1)
            parts = data[0].split()
            self.assertAlmostEqual(float(parts[0]), 400.0 / 3600.0, places=5)
            self.assertAlmostEqual(float(parts[1]), 600.0 / 3600.0, places=5)
            self.assertAlmostEqual(float(parts[2]), N26.field_ra_deg, places=4)
            self.assertAlmostEqual(float(parts[3]), N26.field_dec_deg, places=5)
            self.assertAlmostEqual(float(parts[4]), N26.epoch_jd[0], places=4)
            self.assertEqual(parts[6], "HST.csv")
            self.assertEqual(parts[7], "N26.eff")

    def test_keplerian_plant_on_n26_field(self):
        path = N26_ROOT / "characterization" / "HST.csv"
        obs = grid_bias.parse_jpl_horizons_icrf(path, N26.epoch_jd[0])
        got = grid_bias.keplerian_at_radec_r(
            44.0, 0.05, 3.0, N26.field_ra_deg, N26.field_dec_deg,
            43.5, obs, f_sign=1.0, node_index=0,
        )
        self.assertIsNotNone(got)
        a, e, inc, node, peri, M = got
        ra, dec = grid_bias.apparent_radec_deg(a, e, inc, node, peri, M, obs)
        sep = grid_bias.sky_separation_deg(ra, dec, N26.field_ra_deg, N26.field_dec_deg)
        self.assertLess(sep * 60.0, 0.1)

    def test_aimed_samples_fill_the_rectangle_not_a_square(self):
        path = N26_ROOT / "characterization" / "HST.csv"
        obs = grid_bias.parse_jpl_horizons_icrf(path, N26.epoch_jd[0])
        rng = np.random.default_rng(9)
        half_w = 0.5 * N26.mosaic_width_deg
        half_h = 0.5 * N26.mosaic_height_deg
        ras, decs = [], []
        for _ in range(30):
            el = grid_bias.sample_aimed_elements(43.5, 0.04, 2.5, obs, rng)
            self.assertIsNotNone(el)
            ra, dec = grid_bias.apparent_radec_deg(43.5, 0.04, *el, obs)
            self.assertLessEqual(abs(ra - N26.field_ra_deg), half_w + 1e-3)
            self.assertLessEqual(abs(dec - N26.field_dec_deg), half_h + 1e-3)
            ras.append(ra)
            decs.append(dec)
        self.assertGreater(max(decs) - min(decs), max(ras) - min(ras))


class NapierCatalog(unittest.TestCase):
    def test_four_objects_not_fv53_bg91(self):
        path = N26_ROOT / "data" / "n26_detections.csv"
        with path.open() as fh:
            rows = list(csv.DictReader(fh))
        names = [r["name"] for r in rows]
        self.assertEqual(names, ["2003BF91", "2003BH91", "2003ABCD", "2003WXYZ"])
        self.assertNotIn("2003BG91", names)
        self.assertNotIn("2000FV53", names)
        for row in rows:
            self.assertGreater(float(row["m_stmag"]), 27.5)
            self.assertLess(float(row["m_stmag"]), 29.5)
            self.assertGreater(float(row["d_bary"]), 40.0)
            self.assertLess(float(row["i"]), 5.0)

    def test_load_detections_uses_catalog_ifree_and_stmag(self):
        rows = grid_bias.load_detections(
            N26_ROOT / "data" / "n26_detections.csv"
        )
        self.assertEqual(len(rows), 4)
        by_name = {r["name"]: r for r in rows}
        self.assertAlmostEqual(by_name["2003BF91"]["ifree"], 1.35)
        self.assertAlmostEqual(by_name["2003ABCD"]["a"], 47.5)
        self.assertLess(by_name["2003BF91"]["Hx"], 12.0)
        self.assertGreater(by_name["2003WXYZ"]["Hx"], by_name["2003BF91"]["Hx"])

    def test_n26_script_binds_survey_defaults(self):
        self.assertAlmostEqual(grid_bias.FIELD_RA_DEG, N26.field_ra_deg)
        self.assertAlmostEqual(grid_bias.H_COLOR_OFFSET, -0.3)
        self.assertEqual(grid_bias.EPOCH_JD, N26.epoch_jd)


if __name__ == "__main__":
    unittest.main()
