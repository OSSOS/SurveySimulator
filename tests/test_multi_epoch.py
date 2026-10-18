"""Multi-epoch survey root discovery and detections_required."""
from __future__ import annotations

import pathlib
import shutil
import tempfile
import unittest

from ossssim import SurveyCharacterization, read_survey_conf

script_directory = pathlib.Path(__file__).parent.resolve()
CFEPS = script_directory / 'data' / 'Surveys' / 'CFEPS'
MINI = script_directory.parent / 'docs' / 'examples' / 'mini_survey'


def _copy_epoch(src: pathlib.Path, dest: pathlib.Path) -> None:
    dest.mkdir(parents=True, exist_ok=True)
    for name in ('pointings.list', '2013AE.eff', '2013AO.eff'):
        shutil.copy2(src / name, dest / name)


class SingleEpochRootTest(unittest.TestCase):
    def test_cfeps_flat(self):
        if not CFEPS.is_dir():
            raise unittest.SkipTest(f'missing {CFEPS}')
        s = SurveyCharacterization.from_directory(str(CFEPS))
        self.assertEqual(s.n_epochs, 1)
        self.assertEqual(s.detections_required, 1)
        self.assertEqual(s.by_index[0].survey, 'CFEPS')
        self.assertTrue(s.by_index[0].key.startswith('CFEPS/'))


class MultiEpochRootTest(unittest.TestCase):
    def setUp(self):
        if not MINI.is_dir():
            raise unittest.SkipTest(f'missing {MINI}')
        self.tmp = tempfile.TemporaryDirectory()
        self.root = pathlib.Path(self.tmp.name) / 'MySurvey'
        self.root.mkdir()
        _copy_epoch(MINI, self.root / 'visitOne')
        _copy_epoch(MINI, self.root / '20260610')

    def tearDown(self):
        self.tmp.cleanup()

    def test_discovers_children_and_names_from_root(self):
        s = SurveyCharacterization.from_directory(str(self.root))
        self.assertEqual(s.n_epochs, 2)
        self.assertEqual(s.detections_required, 1)
        self.assertGreater(len(s), 0)
        for p in s:
            self.assertEqual(p.survey, 'MySurvey')
            self.assertTrue(p.key.startswith('MySurvey/'))
            self.assertNotIn('visitOne', p.key)
            self.assertNotIn('20260610', p.key)

    def test_survey_conf_detections_required(self):
        (self.root / 'survey.conf').write_text(
            '# multi-epoch AND\ndetections_required = 2\n'
        )
        self.assertEqual(read_survey_conf(self.root), 2)
        s = SurveyCharacterization.from_directory(str(self.root))
        self.assertEqual(s.n_epochs, 2)
        self.assertEqual(s.detections_required, 2)

    def test_detections_required_out_of_range_fails(self):
        (self.root / 'survey.conf').write_text('detections_required = 3\n')
        with self.assertRaises(RuntimeError):
            SurveyCharacterization.from_directory(str(self.root))

    def test_empty_root_without_children_fails(self):
        empty = pathlib.Path(self.tmp.name) / 'EmptySur'
        empty.mkdir()
        with self.assertRaises(FileNotFoundError):
            SurveyCharacterization.from_directory(str(empty))

    def test_root_pointings_list_wins_over_children(self):
        # Put pointings.list at root as well — single-epoch mode.
        shutil.copy2(MINI / 'pointings.list', self.root / 'pointings.list')
        for name in ('2013AE.eff', '2013AO.eff'):
            shutil.copy2(MINI / name, self.root / name)
        s = SurveyCharacterization.from_directory(str(self.root))
        self.assertEqual(s.n_epochs, 1)
        self.assertEqual(s.by_index[0].survey, 'MySurvey')


if __name__ == '__main__':
    unittest.main()
