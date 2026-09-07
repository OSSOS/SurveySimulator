"""Tests for SurveyCharacterization area / fill / depth coverage API."""
import pathlib
import unittest

import numpy as np

from ossssim import Characterizations, SurveyCharacterization

script_directory = pathlib.Path(__file__).parent.resolve()
CFEPS = Characterizations.surveys['CFEPS']


class SurveyCharacterizationTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.survey = SurveyCharacterization.from_directory(CFEPS)

    def test_loads_pointings(self):
        self.assertGreater(len(self.survey), 0)
        self.assertIn('L3f-smooth', self.survey.keys())
        self.assertIn('L3h-smooth', self.survey.pointings)

    def test_duplicate_eff_keys_disambiguated(self):
        # L3w has two pointings sharing L3w-smooth.eff
        keys = self.survey.keys()
        self.assertTrue(
            any(k.startswith('L3w-smooth') for k in keys),
            msg=f'expected L3w-smooth keys among {keys[:10]}...',
        )
        l3w = [k for k in keys if k.startswith('L3w-smooth')]
        self.assertGreaterEqual(len(l3w), 2)
        self.assertIn('L3w-smooth#0', l3w)
        self.assertIn('L3w-smooth#1', l3w)

    def test_area_and_fill(self):
        p = self.survey['L3f-smooth']
        self.assertGreater(p.area_deg2, 0.0)
        self.assertGreater(p.fill_factor, 0.0)
        self.assertLessEqual(p.fill_factor, 1.0)
        self.assertAlmostEqual(p.fill_factor, 0.80, places=2)

    def test_efficiency_and_effective_area(self):
        p = self.survey['L3h-smooth']
        eta = p.efficiency(23.0)
        self.assertGreater(eta, 0.0)
        self.assertLessEqual(eta, 1.0)
        a_eff = p.effective_area(23.0)
        self.assertAlmostEqual(
            a_eff, p.area_deg2 * p.fill_factor * eta, places=10
        )

    def test_coverage_vs_magnitude_decreases(self):
        mags = np.arange(21.0, 26.01, 0.5)
        m, area = self.survey.coverage_vs_magnitude(mags)
        self.assertEqual(list(m), list(mags))
        self.assertTrue(np.all(area >= 0.0))
        # Fainter end should not be larger than bright end for survey coverage
        self.assertLessEqual(area[-1], area[0] + 1e-9)
        # Check consistency with single-pointing sum at one mag
        mid = 23.0
        _, a_mid = self.survey.coverage_vs_magnitude([mid])
        expected = sum(p.effective_area(mid) for p in self.survey)
        self.assertAlmostEqual(float(a_mid[0]), float(expected), places=8)

    def test_subset_pointing_ids(self):
        mags = np.array([22.0, 24.0])
        _, a_all = self.survey.coverage_vs_magnitude(mags)
        _, a_one = self.survey.coverage_vs_magnitude(
            mags, pointing_ids=['L3f-smooth']
        )
        self.assertTrue(np.all(a_one <= a_all + 1e-9))
        p = self.survey['L3f-smooth']
        np.testing.assert_allclose(
            a_one, p.effective_area(mags), rtol=1e-10
        )


if __name__ == '__main__':
    unittest.main()
