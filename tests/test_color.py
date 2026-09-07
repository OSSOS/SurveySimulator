import unittest
from ossssim import ModelFileOld
from astropy import units as u
import numpy
import pathlib

script_directory = pathlib.Path(__file__).parent.resolve()

class TestOldStyleColorMapping(unittest.TestCase):

    def setUp(self):
        self.old_stype_model_file = f"{script_directory}/data/test_model.dat"
        self.old_style_color_list = numpy.array([0.0, -0.7, -1.2, -1.7, 0.8, 0.5, 0.1, -0.8, -1.2, 0.0]) * u.mag
        self.model_file = ModelFileOld(filename=self.old_stype_model_file)
        self.old_style_color_header = "# Colors = 0.0d0 -0.70d0 -1.2d0 -1.7d0 0.8d0 0.5d0  0.1d0 -0.8d0 -1.2d0 0.0d0"

    def test_color_init(self):
        color_dict = self.model_file.colors('default', 'g')
        self.assertAlmostEqual(color_dict['g-g'], self.old_style_color_list[0])
        self.assertAlmostEqual(color_dict['r-g'], self.old_style_color_list[1])

    def test_color_transform(self):
        color_dict = self.model_file.colors('default', 'r')
        self.assertAlmostEqual(color_dict['g-r'], self.old_style_color_list[0]-self.old_style_color_list[1])
        self.assertAlmostEqual(color_dict['r-r'], self.old_style_color_list[0])

    def tearDown(self):
        self.model_file.close()

    def test_color_list(self):
        color_list = self.model_file.colors.colors_list('default', 'g') * u.mag
        self.assertEqual(self.model_file.colors('default', 'g')['g-g'], self.old_style_color_list[0])
        self.assertEqual(self.model_file.colors('default', 'g')['g-g'], color_list[ord('g')-ord('A')+1])
        self.assertEqual(self.old_style_color_list[0], color_list[ord('g')-ord('A')+1])
        self.assertEqual(self.old_style_color_list[1], color_list[ord('r')-ord('A')+1])


if __name__ == '__main__':
    unittest.main()