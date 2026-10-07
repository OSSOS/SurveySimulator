"""Tests for streaming ECSV ModelFile reads (no upfront full-table load)."""
import os
import unittest
from tempfile import NamedTemporaryFile

from astropy import units
from astropy.time import Time

from ossssim import ModelFile, ModelOutputFile, PhotSpec, definitions
from ossssim.models import ModelFileEcsv


def _sample_row(a_au=40.0, comp='Ring'):
    return {
        'a': a_au * units.au,
        'e': 0.1,
        'inc': 10 * units.deg,
        'node': 20 * units.deg,
        'peri': 30 * units.deg,
        'M': 40 * units.deg,
        'H': 8 * units.mag,
        'q': 36 * units.au,
        'comp': comp,
        'j': 0,
        'k': 0,
        'x': 1 * units.au,
        'y': 2 * units.au,
        'z': 3 * units.au,
        'flag': 1,
    }


def _write_minimal_ecsv(path, datatype_lines, colnames, data_lines, meta_lines=None):
    """Write a tiny hand-crafted ECSV for streaming edge-case tests."""
    meta_lines = meta_lines or [
        '- {Seed: 1}',
        '- {Model_Band: r}',
    ]
    with open(path, 'w') as f_obj:
        f_obj.write('# %ECSV 1.0\n')
        f_obj.write('# ---\n')
        f_obj.write('# datatype:\n')
        for line in datatype_lines:
            f_obj.write(f'# - {line}\n')
        f_obj.write("# delimiter: ','\n")
        f_obj.write('# meta: !!omap\n')
        for line in meta_lines:
            f_obj.write(f'# {line}\n')
        f_obj.write('# schema: astropy-2.0\n')
        f_obj.write(','.join(colnames) + '\n')
        for line in data_lines:
            f_obj.write(line + '\n')


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


class ModelFileEcsvCodexRegressionTest(unittest.TestCase):
    """Regression coverage for Codex review findings on PR #11."""

    def setUp(self):
        self._tmp = NamedTemporaryFile(suffix='.ecsv', delete=False)
        self._tmp.close()
        self.filename = self._tmp.name
        if os.path.exists(self.filename):
            os.remove(self.filename)

    def tearDown(self):
        if os.path.exists(self.filename):
            os.remove(self.filename)

    def test_randomize_one_row_file_does_not_hang(self):
        writer = ModelOutputFile(
            filename=self.filename,
            seed=7,
            epoch=definitions.Neptune['Epoch'],
            longitude_neptune=definitions.Neptune['Longitude'],
            colors=PhotSpec(),
            model_band='r',
        )
        writer.write_row(_sample_row(a_au=42.5))

        model = ModelFile(self.filename, randomize=True)
        for _ in range(5):
            row = next(model)
            self.assertAlmostEqual(row['a'].to('au').value, 42.5, places=4)
            self.assertEqual(row['comp'], 'Ring')
        self.assertIsNone(model._table)
        model.close()

    def test_stream_honors_ecsv_schema_units(self):
        _write_minimal_ecsv(
            self.filename,
            datatype_lines=[
                "{name: a, unit: AU, datatype: float64}",
                "{name: e, datatype: float64}",
                "{name: inc, unit: rad, datatype: float64}",
                "{name: comp, datatype: string}",
            ],
            colnames=['a', 'e', 'inc', 'comp'],
            data_lines=['40.0,0.1,1.0,Ring'],
            meta_lines=[
                '- {Seed: 1}',
                '- {Model_Band: r}',
                f"- Epoch: {definitions.Neptune['Epoch'].jd}",
            ],
        )
        model = ModelFile(self.filename)
        row = next(model)
        self.assertIsNone(model._table)
        self.assertEqual(row['inc'].unit, units.rad)
        self.assertAlmostEqual(row['inc'].value, 1.0, places=6)
        self.assertAlmostEqual(row['a'].to('au').value, 40.0, places=6)
        model.close()

    def test_stream_parses_quoted_delimiter_fields(self):
        _write_minimal_ecsv(
            self.filename,
            datatype_lines=[
                "{name: a, unit: AU, datatype: float64}",
                "{name: e, datatype: float64}",
                "{name: comp, datatype: string}",
                "{name: j, datatype: int32}",
                "{name: k, datatype: int32}",
            ],
            colnames=['a', 'e', 'comp', 'j', 'k'],
            data_lines=['40.0,0.1,"A,B",5,1'],
            meta_lines=[
                '- {Seed: 1}',
                '- {Model_Band: r}',
            ],
        )
        model = ModelFile(self.filename)
        row = next(model)
        self.assertIsNone(model._table)
        self.assertEqual(row['comp'], 'A,B')
        self.assertEqual(row['j'], 5)
        self.assertEqual(row['k'], 1)
        self.assertAlmostEqual(row['a'].to('au').value, 40.0, places=6)
        model.close()


if __name__ == '__main__':
    unittest.main()
