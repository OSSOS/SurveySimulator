import unittest

import numpy
from astropy.units import Quantity

import osssim
from osssim import definitions

class OSSSIMTest(unittest.TestCase):
    def setUp(self):
        self.model = osssim.SSimModelFile('data/test_model.dat')
        self.detect_filename = 'test_detect.dat'
        result = osssim.SSimModelFile('data/test_detect.dat')
        self.result_row = next(iter(result))
        self.seed = int(result.header['Seed'])
        self.osssim = osssim.OSSSSim('data/CFEPS')

    def test_simulate(self):
        # loop over the model file until we have a detection and then compare the detected row values to the test row
        for row in self.model:
            row = self.osssim.simulate(row, seed=self.seed, epoch=self.model.epoch)
            if row['flag'] > 0:
                for key in self.result_row:
                    test_value = self.result_row[key]
                    result_value = row[key]
                    if isinstance(test_value, Quantity):
                        test_value = test_value.to(osssim.definitions.colunits[key]).value
                        result_value = result_value.to(osssim.definitions.colunits[key]).value
                    self.assertAlmostEqual(test_value, result_value, 4)
                break

if __name__ == '__main__':
    unittest.main()
