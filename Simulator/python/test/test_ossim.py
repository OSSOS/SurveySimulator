import unittest
from ssim_util import SSimModelFile
from ossim import OSSSIM

class OSSIMTest(unittest.TestCase):
    def setUp(self):
        self.model = SSimModelFile('data/L7model-3.0-9.0')
        self.detect = SSimModelFile('data/ReadModelFromFile-check.dat')
        self.seed = 123456789
        self.ossim = OSSSIM(self.model, 'data/CFEPS')

    def test_simulate(self):
        self.ossim.simulate('data/simulate_test_out.txt', seed=self.seed, ntrack=1)

if __name__ == '__main__':
    unittest.main()