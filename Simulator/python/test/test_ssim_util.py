import unittest
from ssim_util import SSimModelFile
from astropy import units

class SSimModelFileTest(unittest.TestCase):

    def setUp(self):
        self.epoch = 2453157.50000
        self.lambda_neptune = 5.489*units.radian
        self.colors = [ float(x.replace("d","e")) for x in "0.0d0 -0.70d0 -1.2d0 -1.7d0 0.8d0 0.5d0  0.1d0 -0.8d0 -1.2d0 0.0d0".split()]
        self.colnames = "  a      e     inc     node    peri    m       H      dist     comp      j  k".split()
        self.row = {'a':     88.253*units.au,
                    'e':      0.564,
                    'inc':    8.807*units.degree,
                    'node': 316.126*units.degree,
                    'peri': 286.639*units.degree,
                    'm':     60.537*units.degree,
                    'H':      8.05*units.mag,
                    'dist':  90.70*units.au,
                    'comp': 'resonant',
                    'j':     5,
                    'k':     1}
        self.model = SSimModelFile('data/L7model-3.0-9.0')


    def test_epoch(self):
        self.assertAlmostEqual(self.epoch, self.model.epoch.jd)

    def test_lambda_neptune(self):
        self.assertAlmostEqual(self.lambda_neptune, self.model.lambda_neptune)

    def test_colors(self):
        self.assertEqual(self.colors, self.model.colors)

    def test_colnames(self):
        self.assertEqual(self.colnames, self.model.colnames)

    def test_row(self):
        for row in self.model:
            for key in self.row:
                self.assertEqual(self.row[key], row[key])
            break

if __name__ == '__main__':
    unittest.main()