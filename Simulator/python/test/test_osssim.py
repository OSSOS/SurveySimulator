import unittest
import osssim

class OSSIMTest(unittest.TestCase):
    def setUp(self):
        self.model = osssim.SSimModelFile('data/test_model.dat')
        result = osssim.SSimModelFile('data/test_detect.dat')
        self.result_row = next(iter(result))
        self.seed = int(result.header['Seed'])
        self.osssim = osssim.OSSSIM(self.model, 'data/CFEPS')

    def test_simulate(self):
        self.osssim.simulate('simulate_test_out.txt', seed=self.seed, ntrack=0)
        data = next(iter(osssim.SSimModelFile('simulate_test_out.txt')))
        for key in self.result_row:
            self.assertEqual(self.result_row[key], data[key])

if __name__ == '__main__':
    unittest.main()
    