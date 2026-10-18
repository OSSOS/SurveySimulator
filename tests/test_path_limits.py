"""Tests for path_len / name_len / block_len and survey/block detection keys."""
from __future__ import annotations

import os
import pathlib
import shutil
import tempfile
import unittest

from ossssim import SurveyCharacterization, OSSSSim, ModelFile

script_directory = pathlib.Path(__file__).parent.resolve()
CFEPS = script_directory / 'data' / 'Surveys' / 'CFEPS'
MINI = script_directory.parent / 'docs' / 'examples' / 'mini_survey'


def _link_tree(src: pathlib.Path, dst: pathlib.Path) -> None:
    """Hard-link or copy survey files into dst (pointings.list + .eff)."""
    dst.mkdir(parents=True, exist_ok=True)
    for item in src.iterdir():
        if item.is_file() and (
            item.name == 'pointings.list'
            or item.suffix.lower() == '.eff'
            or item.suffix.lower() == '.csv'
        ):
            target = dst / item.name
            try:
                os.link(item, target)
            except OSError:
                shutil.copy2(item, target)


class SurveyBlockKeyTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not CFEPS.is_dir():
            raise unittest.SkipTest(f'CFEPS fixture missing: {CFEPS}')
        cls.survey = SurveyCharacterization.from_directory(str(CFEPS))

    def test_detection_key_format(self):
        p = self.survey['L3f-smooth']
        self.assertEqual(p.survey, 'CFEPS')
        self.assertEqual(p.block, 'L3f-smooth')
        self.assertEqual(p.key, 'CFEPS/L3f-smooth')
        self.assertTrue(p.efnam.endswith('.eff'))

    def test_simulate_writes_survey_block_key(self):
        model = ModelFile(f'{script_directory}/data/test_model.dat')
        detect = ModelFile(f'{script_directory}/data/test_detect.dat')
        seed = int(detect.header['Seed'][0])
        detect.close()
        sim = OSSSSim(characterization_directory=str(CFEPS), seed=seed)
        found = None
        for row in model:
            result = sim.simulate(
                row, epoch=model.epoch, colors=model.colors, model_band=model.model_band
            )
            if result['flag'] > 0:
                found = result
                break
        model.close()
        self.assertIsNotNone(found)
        survey_key = found['Survey']
        if isinstance(survey_key, bytes):
            survey_key = survey_key.decode('utf-8')
        survey_key = str(survey_key).strip()
        self.assertIn('/', survey_key)
        survey, block = survey_key.split('/', 1)
        self.assertEqual(survey, 'CFEPS')
        self.assertTrue(len(block) > 0)
        self.assertLessEqual(len(survey_key), 65)


class PathAndNameLimitTest(unittest.TestCase):
    def test_long_characterization_path_loads(self):
        if not MINI.is_dir():
            raise unittest.SkipTest(f'mini_survey missing: {MINI}')
        # Nested path well above the old 100-char limit, under path_len=2048.
        with tempfile.TemporaryDirectory() as tmp:
            deep = pathlib.Path(tmp)
            for part in [f'seg{i:02d}_xxxxxxxx' for i in range(12)]:
                deep = deep / part
            deep = deep / 'mini_survey'
            _link_tree(MINI, deep)
            self.assertGreater(len(str(deep)), 100)
            survey = SurveyCharacterization.from_directory(str(deep))
            self.assertGreater(len(survey), 0)
            self.assertEqual(survey.by_index[0].survey, 'mini_survey')
            self.assertIn('/', survey.by_index[0].key)

    def test_survey_name_too_long_fails(self):
        if not MINI.is_dir():
            raise unittest.SkipTest(f'mini_survey missing: {MINI}')
        with tempfile.TemporaryDirectory() as tmp:
            # name_len = 32; 33-char basename must fail
            long_name = 'S' * 33
            dest = pathlib.Path(tmp) / long_name
            _link_tree(MINI, dest)
            with self.assertRaises(RuntimeError):
                SurveyCharacterization.from_directory(str(dest))

    def test_missing_eff_fails_hard(self):
        if not MINI.is_dir():
            raise unittest.SkipTest(f'mini_survey missing: {MINI}')
        with tempfile.TemporaryDirectory() as tmp:
            dest = pathlib.Path(tmp) / 'BrokenSur'
            _link_tree(MINI, dest)
            # Remove one referenced .eff so load must fail (no skip-and-continue).
            effs = list(dest.glob('*.eff'))
            self.assertTrue(effs)
            effs[0].unlink()
            with self.assertRaises(RuntimeError):
                SurveyCharacterization.from_directory(str(dest))


if __name__ == '__main__':
    unittest.main()
