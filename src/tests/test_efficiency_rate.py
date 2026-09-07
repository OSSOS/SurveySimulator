"""Tests for rate-dependent m0 in single hyperbolic-tangent efficiency."""
import math
import unittest

from ossssim.characterization import SingleTanhParam


class SingleTanhRateDependenceTest(unittest.TestCase):
    def setUp(self):
        self.eta0 = 0.90
        self.m00 = 24.00
        self.w = 0.50
        self.alpha = 0.10
        self.mag = 24.00
        self.param = SingleTanhParam(
            A=self.eta0, M_0=self.m00, w=self.w, alpha_rate=self.alpha
        )

    @staticmethod
    def expected(eta0, m00, w, alpha, mag, rate_asphr):
        m0 = m00 + alpha * rate_asphr
        return eta0 / 2.0 * (1.0 - math.tanh((mag - m0) / w))

    def test_rate_raises_m0_and_efficiency(self):
        eta_slow = self.param.efficiency(self.mag, rate_asphr=1.0)
        eta_fast = self.param.efficiency(self.mag, rate_asphr=3.0)
        self.assertAlmostEqual(
            eta_slow,
            self.expected(self.eta0, self.m00, self.w, self.alpha, self.mag, 1.0),
            places=12,
        )
        self.assertAlmostEqual(
            eta_fast,
            self.expected(self.eta0, self.m00, self.w, self.alpha, self.mag, 3.0),
            places=12,
        )
        # Faster motion -> brighter (larger) m0 -> higher efficiency at fixed mag
        self.assertGreater(eta_fast, eta_slow)

    def test_zero_alpha_is_rate_independent(self):
        legacy = SingleTanhParam(A=self.eta0, M_0=self.m00, w=self.w, alpha_rate=0.0)
        eta1 = legacy.efficiency(self.mag, rate_asphr=1.0)
        eta3 = legacy.efficiency(self.mag, rate_asphr=3.0)
        expected = self.expected(self.eta0, self.m00, self.w, 0.0, self.mag, 0.0)
        self.assertAlmostEqual(eta1, expected, places=12)
        self.assertAlmostEqual(eta3, expected, places=12)

    def test_round_trip_string_omits_zero_alpha(self):
        legacy = SingleTanhParam(A=0.9, M_0=24.0, w=0.5, alpha_rate=0.0)
        self.assertEqual(str(legacy), "single_param= 0.90 24.00 0.50")
        with_alpha = SingleTanhParam(A=0.9, M_0=24.0, w=0.5, alpha_rate=0.1)
        self.assertIn("0.1", str(with_alpha))
        parsed = SingleTanhParam.from_string(str(with_alpha))
        self.assertAlmostEqual(parsed.alpha_rate, 0.1, places=6)

    def test_from_string_legacy_three_params(self):
        parsed = SingleTanhParam.from_string("single_param= 0.90 24.00 0.50")
        self.assertAlmostEqual(parsed.A, 0.9)
        self.assertAlmostEqual(parsed.M_0, 24.0)
        self.assertAlmostEqual(parsed.w, 0.5)
        self.assertAlmostEqual(parsed.alpha_rate, 0.0)


if __name__ == "__main__":
    unittest.main()
