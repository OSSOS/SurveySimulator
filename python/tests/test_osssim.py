import unittest
from astropy.units import Quantity
import ossssim
from tempfile import NamedTemporaryFile

class OSSSIMTest(unittest.TestCase):
    def setUp(self):
        self.model = ossssim.ModelFile('data/test_model.dat')
        self.detect_fobj = NamedTemporaryFile()
        self.detect_filename = self.detect_fobj.name
        result = ossssim.ModelFile('data/test_detect.dat')
        self.result_row = next(iter(result))
        self.seed = int(result.header['Seed'])
        self.osssim = ossssim.OSSSSim('data/Surveys/CFEPS')

    def test_simulate(self):
        # loop over the model file until we have a detection and then compare the detected row values to the test row
        for row in self.model:
            result_row = self.osssim.simulate(row, seed=self.seed, epoch=self.model.epoch)
            if result_row['flag'] > 0:
                for key in row:
                    test_value = row[key]
                    result_value = result_row[key]
                    if isinstance(test_value, Quantity):
                        test_value = test_value.to(ossssim.definitions.colunits[key]).value
                        result_value = result_value.to(ossssim.definitions.colunits[key]).value
                    self.assertAlmostEqual(test_value, result_value, 4)
                break


if __name__ == '__main__':
    unittest.main()
