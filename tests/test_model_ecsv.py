"""Tests for streaming ECSV ModelFile reads (no upfront full-table load)."""
import os
import unittest
from tempfile import NamedTemporaryFile

from astropy import units
from astropy.time import Time

from ossssim import ModelFile, ModelOutputFile, PhotSpec, definitions
from ossssim.models import ModelFileEcsv


def _sample_row(a_au=40.0):
    return {
        'a': a_au * units.au,
        'e': 0.1,
        'inc': 10 * units.deg,
        'node': 20 * units.deg,
        'peri': 30 * units.deg,
        'M': 40 * units.deg,
        'H': 8 * units.mag,
        'q': 36 * units.au,
        'comp': 'Ring',
        'j': 0,
        'k': 0,
        'x': 1 * units.au,
        'y': 2 * units.au,
        'z': 3 * units.au,
        'flag': 1,
    }


class ModelFileEcsvStreamingTest(unittest.TestCase):

    def setUp(self):
        self._tmp = NamedTemporaryFile(suffix='.ecsv', delete=False)
        self._tmp.close()
        self.filename = self._tmp.name
        if os.path.exists(self.filename):
            os.remove(self.filename)

        self.epoch = definitions.Neptune['Epoch']
        self.seed = 424242
        self.model_band = 'r'
        self.colors = PhotSpec()
        writer = ModelOutputFile(
            filename=self.filename,
            seed=self.seed,
            epoch=self.epoch,
            longitude_neptune=definitions.Neptune['Longitude'],
            colors=self.colors,
            model_band=self.model_band,
        )
        self.n_rows = 8
        for i in range(self.n_rows):
            writer.write_row(_sample_row(a_au=40.0 + i))

    def tearDown(self):
        if os.path.exists(self.filename):
            os.remove(self.filename)

    def test_factory_selects_ecsv(self):
        model = ModelFile(self.filename)
        self.assertIsInstance(model, ModelFileEcsv)
        model.close()

    def test_header_without_loading_table(self):
        model = ModelFile(self.filename)
        self.assertIsNone(model._table)
        self.assertEqual(model.seed, self.seed)
        self.assertEqual(model.model_band, self.model_band)
        self.assertIsInstance(model.epoch, Time)
        self.assertAlmostEqual(model.epoch.jd, self.epoch.jd)
        self.assertIn('default', model.colors.colors)
        self.assertIsNone(model._table)
        model.close()

    def test_stream_does_not_load_table(self):
        model = ModelFile(self.filename)
        row = next(model)
        self.assertIsNone(model._table)
        self.assertEqual(row['comp'], 'Ring')
        self.assertAlmostEqual(row['a'].to('au').value, 40.0, places=5)
        self.assertAlmostEqual(row['e'], 0.1, places=5)
        model.close()

    def test_streamed_rows_match_table(self):
        streamed = ModelFile(self.filename)
        rows = [row for row in streamed]
        self.assertIsNone(streamed._table)
        streamed.close()

        full = ModelFile(self.filename)
        table = full.table
        self.assertEqual(len(rows), len(table))
        for i, row in enumerate(rows):
            self.assertAlmostEqual(row['a'].to('au').value,
                                   table['a'][i].to('au').value, places=4)
            self.assertEqual(row['comp'], table['comp'][i])
            self.assertEqual(row['j'], table['j'][i])
        full.close()

    def test_randomize_returns_valid_rows_without_loading_table(self):
        model = ModelFile(self.filename, randomize=True)
        self.assertTrue(model.randomize)
        colnames = set(model.colnames)
        for _ in range(10):
            row = next(model)
            self.assertEqual(set(row.keys()), colnames)
            self.assertEqual(row['comp'], 'Ring')
            self.assertIsNotNone(row['a'])
            a_val = row['a'].to('au').value
            self.assertGreaterEqual(a_val, 40.0)
            self.assertLess(a_val, 40.0 + self.n_rows)
        self.assertIsNone(model._table)
        model.close()


if __name__ == '__main__':
    unittest.main()
