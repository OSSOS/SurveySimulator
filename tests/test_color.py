"""
Unit tests for PhotSpec color / bandpass indexing and H-magnitude band transfer.

ASCII color array layout (shared with Fortran filter_to_index):
  Python index i = ord(L) - ord('A')
  Fortran index  = IACHAR(L) - IACHAR('A') + 1
  so Python[i] maps to Fortran(i+1) across f2py.

Photometric contract (Detos1 / Driver):
  - Input H (hx) is in the model bandpass ``x``
  - color(filter) = mag_filter - mag_x
  - Detection uses H_filter = H_x + color(filter)
  - Returned m_int is converted back to band ``x``
  - Returned m_rand / h_rand remain in the discovery filter
  - Column H itself is never rewritten (stays model-band)
"""
from __future__ import annotations

import pathlib
import unittest

import numpy
from astropy import units as u

from ossssim.color import PhotSpec
from ossssim import ModelFileOld

script_directory = pathlib.Path(__file__).parent.resolve()

# Classic optical + wide + near-IR letters used by defaults / surveys
COMMON_BANDS = ['g', 'r', 'i', 'z', 'u', 'V', 'B', 'R', 'I', 'w', 'J', 'H']


def py_index(band: str) -> int:
    return ord(band) - ord('A')


def fortran_index(band: str) -> int:
    """1-based index matching F95/ioutils.f95 filter_to_index."""
    return ord(band) - ord('A') + 1


class TestAsciiBandpassIndexing(unittest.TestCase):
    """Single-letter bandpass → color-array slot mapping."""

    def test_array_length_covers_A_through_z(self):
        self.assertEqual(ord('z') - ord('A') + 1, 58)
        colors = PhotSpec().colors_list('default', 'g')
        self.assertEqual(len(colors), 58)

    def test_common_band_slots(self):
        # Spot-check a few known ASCII offsets
        self.assertEqual(py_index('A'), 0)
        self.assertEqual(py_index('g'), 38)
        self.assertEqual(py_index('r'), 49)
        self.assertEqual(py_index('z'), 57)
        self.assertEqual(fortran_index('g'), 39)
        self.assertEqual(fortran_index('r'), 50)

    def test_colors_list_stores_at_python_index(self):
        phot = PhotSpec()
        color_list = phot.colors_list('default', 'g')
        self.assertAlmostEqual(color_list[py_index('g')], 0.0)
        self.assertAlmostEqual(color_list[py_index('r')], -0.7)
        self.assertAlmostEqual(color_list[py_index('i')], -1.2)
        self.assertAlmostEqual(color_list[py_index('w')], -0.8)
        self.assertAlmostEqual(color_list[py_index('J')], -1.5)
        self.assertAlmostEqual(color_list[py_index('H')], -1.6)
        # Neighbour of r must not accidentally hold the r-g value (guards +1 bugs)
        self.assertAlmostEqual(color_list[py_index('r') + 1], 0.0)

    def test_python_fortran_index_bridge(self):
        for band in COMMON_BANDS:
            self.assertEqual(fortran_index(band), py_index(band) + 1)
            # Fortran index_to_filter(idx) = CHAR(idx + ICHAR('A') - 1)
            self.assertEqual(chr(fortran_index(band) + ord('A') - 1), band)

    def test_reject_out_of_range_bandpass(self):
        bad = PhotSpec(colors={'default': {'$-g': 0.1 * u.mag, 'g-g': 0.0 * u.mag}})
        with self.assertRaises(ValueError):
            bad.colors_list('default', 'g')


class TestModelBandTransform(unittest.TestCase):
    """Rebase color dictionary onto the model H bandpass."""

    def setUp(self):
        self.phot = PhotSpec()

    def test_identity_when_model_band_is_reference(self):
        colors = self.phot.transform_spectral_group_to_model_band('default', 'g')
        self.assertAlmostEqual(colors['g-g'].to(u.mag).value, 0.0)
        self.assertAlmostEqual(colors['r-g'].to(u.mag).value, -0.7)
        self.assertAlmostEqual(colors['i-g'].to(u.mag).value, -1.2)

    def test_rebase_to_r(self):
        # Defaults are g-relative; r-g = -0.7 ⇒ g-r = +0.7, r-r = 0, i-r = -0.5
        colors = self.phot.transform_spectral_group_to_model_band('default', 'r')
        self.assertAlmostEqual(colors['g-r'].to(u.mag).value, 0.7)
        self.assertAlmostEqual(colors['r-r'].to(u.mag).value, 0.0)
        self.assertAlmostEqual(colors['i-r'].to(u.mag).value, -0.5)

    def test_colors_list_after_rebase(self):
        color_list = self.phot.colors_list('default', 'r')
        self.assertAlmostEqual(color_list[py_index('g')], 0.7)
        self.assertAlmostEqual(color_list[py_index('r')], 0.0)
        self.assertAlmostEqual(color_list[py_index('i')], -0.5)


class TestHBandpassTransferMath(unittest.TestCase):
    """
    Document the Detos1 color algebra without calling Fortran.

    color(filter) = mag_filter - mag_x
    H_filter = H_x + color(filter)
    m_int_x = m_int_filter - color(filter)
    """

    def setUp(self):
        self.phot = PhotSpec()
        self.Hx = 8.0  # absolute mag in model band g
        self.colors = self.phot.transform_spectral_group_to_model_band('default', 'g')

    def test_h_converted_to_survey_filter_then_mint_back_to_model(self):
        color_rg = self.colors['r-g'].to(u.mag).value  # -0.7
        H_r = self.Hx + color_rg
        self.assertAlmostEqual(H_r, 7.3)

        # Suppose AppMag returned m_int in r of 22.0; undo color for model-band m_int
        m_int_r = 22.0
        m_int_g = m_int_r - color_rg
        self.assertAlmostEqual(m_int_g, 22.7)

    def test_zero_color_when_survey_matches_model_band(self):
        color_gg = self.colors['g-g'].to(u.mag).value
        self.assertAlmostEqual(color_gg, 0.0)
        self.assertAlmostEqual(self.Hx + color_gg, self.Hx)

    def test_h_column_unchanged_h_rand_in_discovery_band(self):
        """
        Contract check: reported H stays Hx; h_rand is AbsMag(m_rand) in discovery filter.
        With non-zero color, h_rand != Hx even if distances/phase cancel in a toy equal-mag case.
        """
        color_rg = self.colors['r-g'].to(u.mag).value
        Hx = self.Hx
        # Toy: same geometry ⇒ m_rand_r ≈ m_int_r ⇒ h_rand_r ≈ Hx + color_rg
        h_rand_r = Hx + color_rg
        self.assertNotAlmostEqual(h_rand_r, Hx)
        self.assertAlmostEqual(Hx, self.Hx)  # model H never shifted


class TestFromListRoundTrip(unittest.TestCase):

    def test_nonzero_terms_round_trip(self):
        phot = PhotSpec()
        color_list = phot.colors_list('default', 'g')
        rebuilt = PhotSpec.from_list(color_list, 'g')
        original = phot.transform_spectral_group_to_model_band('default', 'g')
        for key, value in original.items():
            if abs(value.to(u.mag).value) < 1e-12:
                continue  # sparse from_list omits exact zeros
            self.assertIn(key, rebuilt.colors['default'])
            self.assertAlmostEqual(
                rebuilt.colors['default'][key].to(u.mag).value,
                value.to(u.mag).value,
            )


class TestOldStyleColorMapping(unittest.TestCase):

    def setUp(self):
        self.old_style_model_file = f"{script_directory}/data/test_model.dat"
        # 9-value header in test_model.dat (no trailing w)
        self.old_style_color_list = numpy.array(
            [0.0, -0.7, -1.2, -1.7, 0.8, 0.5, 0.1, -0.8, -1.2]
        ) * u.mag
        self.model_file = ModelFileOld(filename=self.old_style_model_file)

    def tearDown(self):
        self.model_file.close()

    def test_color_init_from_file(self):
        color_dict = self.model_file.colors('default', 'g')
        self.assertAlmostEqual(color_dict['g-g'], self.old_style_color_list[0])
        self.assertAlmostEqual(color_dict['r-g'], self.old_style_color_list[1])

    def test_color_transform_to_r(self):
        color_dict = self.model_file.colors('default', 'r')
        self.assertAlmostEqual(
            color_dict['g-r'].to(u.mag).value,
            (self.old_style_color_list[0] - self.old_style_color_list[1]).to(u.mag).value,
        )
        self.assertAlmostEqual(color_dict['r-r'].to(u.mag).value, 0.0)

    def test_colors_list_ascii_slots(self):
        color_list = self.model_file.colors.colors_list('default', 'g')
        self.assertAlmostEqual(color_list[py_index('g')], self.old_style_color_list[0].to(u.mag).value)
        self.assertAlmostEqual(color_list[py_index('r')], self.old_style_color_list[1].to(u.mag).value)

    def test_ten_value_old_list_includes_w(self):
        ten = numpy.array([0.0, -0.7, -1.2, -1.7, 0.8, 0.5, 0.1, -0.8, -1.2, -0.8]) * u.mag
        phot = PhotSpec.from_old_style_list(ten)
        colors = phot.transform_spectral_group_to_model_band('default', 'g')
        self.assertAlmostEqual(colors['w-g'].to(u.mag).value, -0.8)


class TestCoreBandColorReporting(unittest.TestCase):
    """Regression for the ic → band / color_offset_array lookup."""

    def test_discovery_band_color_lookup(self):
        phot = PhotSpec()
        color_offset_array = phot.colors_list('default', 'g')
        # Simulate Fortran returning filt_i for 'r'
        ic = fortran_index('r')
        band = chr(ic + ord('A') - 1)
        color = color_offset_array[ord(band) - ord('A')]
        self.assertEqual(band, 'r')
        self.assertAlmostEqual(color, -0.7)
        # The previous off-by-one would read the next slot (0.0)
        self.assertNotAlmostEqual(color_offset_array[ord(band) - ord('A') + 1], -0.7)


class TestFortranFilterIndexAlignment(unittest.TestCase):
    """Optional: compare Python slots to compiled ossssimlib.ioutils."""

    @classmethod
    def setUpClass(cls):
        cls.ossssimlib = None
        try:
            import ossssimlib  # noqa: F401
            cls.ossssimlib = ossssimlib
        except ImportError:
            pass

    def test_filter_to_index_matches_python(self):
        if self.ossssimlib is None:
            self.skipTest('ossssimlib not importable')
        for band in COMMON_BANDS:
            f_idx = int(self.ossssimlib.ioutils.filter_to_index(band))
            self.assertEqual(f_idx, fortran_index(band), msg=band)
            letter = self.ossssimlib.ioutils.index_to_filter(f_idx)
            if isinstance(letter, bytes):
                letter = letter.decode('ascii')
            self.assertEqual(str(letter)[0], band, msg=f'index_to_filter({f_idx})')


if __name__ == '__main__':
    unittest.main()
