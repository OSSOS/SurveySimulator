"""
Survey characterization access: areas, fill factors, and efficiency vs depth.

Example
-------
>>> from ossssim import SurveyCharacterization
>>> import numpy as np
>>> survey = SurveyCharacterization.from_directory('F95/tests/Surveys/CFEPS')
>>> mags = np.arange(21.0, 26.0, 0.1)
>>> m, area = survey.coverage_vs_magnitude(mags)
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple, Union
import numpy as np

import ossssimlib

MagLike = Union[float, np.ndarray, Sequence[float]]


def _as_str(value) -> str:
    if isinstance(value, bytes):
        return value.decode('utf-8', errors='ignore').strip()
    return str(value).strip()


@dataclass
class Pointing:
    """One FoV pointing with geometry, fill factor, and efficiency (block)."""

    id: str
    index: int  # 1-based Fortran index into the loaded survey
    survey: str  # characterization directory basename
    block: str  # .eff stem; unique within survey only
    key: str  # "survey/block" detection key ('/' is a delimiter, not a path)
    efnam: str  # efficiency basename as listed in pointings.list
    ra: float  # radians
    dec: float  # radians
    epoch: float  # JD
    area_deg2: float
    fill_factor: float
    mag_lim: float
    obs_code: int
    rate_mid_asphr: float
    _lib: object = field(repr=False, default=None)

    def efficiency(
        self,
        mag: MagLike,
        rate_asphr: Optional[float] = None,
    ) -> Union[float, np.ndarray]:
        """Evaluate η at magnitude(s). Default rate is the pointing rate_cut midpoint."""
        rate = self.rate_mid_asphr if rate_asphr is None else float(rate_asphr)
        lib = self._lib or ossssimlib.surveysub
        arr = np.asarray(mag, dtype=float)
        if arr.ndim == 0:
            # f90wrap returns intent(out) maglim before the function result
            _maglim, eta = lib.pointing_eta(self.index, float(arr), rate)
            return float(eta)
        etas = np.empty(arr.size, dtype=float)
        lib.pointing_eta_grid(
            self.index, arr.ravel(), arr.size, rate, etas
        )
        return etas.reshape(arr.shape)

    def effective_area(
        self,
        mag: MagLike,
        rate_asphr: Optional[float] = None,
        include_fill: bool = True,
    ) -> Union[float, np.ndarray]:
        """Return area * fill * η(m) (fill omitted when include_fill is False)."""
        eta = self.efficiency(mag, rate_asphr=rate_asphr)
        scale = self.area_deg2 * (self.fill_factor if include_fill else 1.0)
        return scale * eta


def read_survey_conf(directory: Union[str, Path]) -> int:
    """
    Read ``detections_required`` from ``survey.conf`` in the survey root.

    Missing file defaults to 1. Same key=value contract as the Fortran reader.
    """
    path = Path(directory) / 'survey.conf'
    detections_required = 1
    if not path.is_file():
        return detections_required
    for raw in path.read_text().splitlines():
        line = raw.strip()
        if not line or line.startswith('#'):
            continue
        if '=' not in line:
            continue
        key, _, val = line.partition('=')
        if key.strip() == 'detections_required':
            detections_required = int(val.strip().split()[0])
    return detections_required


class SurveyCharacterization:
    """
    Loaded survey characterization root with keyed pointing access.

    Pass the *survey root*. If ``root/pointings.list`` exists it is a
    single-epoch survey; otherwise each immediate child containing
    ``pointings.list`` is an epoch. Survey identity is the root basename.
    Optional ``survey.conf`` sets ``detections_required`` (default 1).

    A *block* is an ``.eff`` stem within that survey. Pointing keys use the
    block name; when several pointings share the same block they are
    disambiguated as ``block#0``, ``block#1``, ... Detection attribution uses
    ``survey/block`` (see ``Pointing.key``). Use ``by_index`` for unambiguous
    iteration.
    """

    def __init__(
        self,
        directory: Path,
        pointings: Dict[str, Pointing],
        by_index: List[Pointing],
        n_epochs: int = 1,
        detections_required: int = 1,
    ):
        self.directory = Path(directory)
        self.pointings = pointings
        self.by_index = by_index
        self.n_epochs = n_epochs
        self.detections_required = detections_required

    @classmethod
    def from_directory(
        cls,
        directory: Union[str, Path],
        lun: int = 21,
    ) -> "SurveyCharacterization":
        """Load a survey root (single- or multi-epoch) via Fortran."""
        directory = Path(directory).resolve()
        if not directory.is_dir():
            raise FileNotFoundError(f"Survey root not found: {directory}")
        has_root_pl = (directory / 'pointings.list').is_file()
        if not has_root_pl:
            kids = [
                p for p in directory.iterdir()
                if p.is_dir() and (p / 'pointings.list').is_file()
            ]
            if not kids:
                raise FileNotFoundError(
                    f"No pointings.list in {directory} or its immediate children"
                )

        lib = ossssimlib.surveysub
        lib.reset_simulator()
        n_sur, ierr = lib.survey_load(str(directory), lun)
        if ierr != 0 or n_sur <= 0:
            raise RuntimeError(
                f"survey_load failed for {directory}: ierr={ierr}, n_sur={n_sur}"
            )
        n_epochs, detections_required = lib.survey_meta()

        # First pass: collect metadata; count duplicate blocks for id disambiguation
        raw = []
        counts: Dict[str, int] = {}
        for i in range(1, n_sur + 1):
            area, fill = lib.pointing_geom(i)
            ra, dec = lib.pointing_center(i)
            survey, block, key, eff_file, epoch, code, mag_lim, rate_mid = (
                lib.pointing_meta(i)
            )
            survey_s = _as_str(survey)
            block_s = _as_str(block)
            key_s = _as_str(key)
            eff_s = _as_str(eff_file)
            counts[block_s] = counts.get(block_s, 0) + 1
            raw.append(
                dict(
                    index=i,
                    survey=survey_s,
                    block=block_s,
                    key=key_s,
                    efnam=eff_s,
                    ra=float(ra),
                    dec=float(dec),
                    epoch=float(epoch),
                    area_deg2=float(area),
                    fill_factor=float(fill),
                    mag_lim=float(mag_lim),
                    obs_code=int(code),
                    rate_mid_asphr=float(rate_mid),
                )
            )

        # Second pass: assign ids (plain block if unique, else block#k)
        seen: Dict[str, int] = {}
        by_index: List[Pointing] = []
        pointings: Dict[str, Pointing] = {}
        for item in raw:
            block = item['block']
            if counts[block] == 1:
                pid = block
            else:
                k = seen.get(block, 0)
                seen[block] = k + 1
                pid = f"{block}#{k}"
            p = Pointing(
                id=pid,
                index=item['index'],
                survey=item['survey'],
                block=item['block'],
                key=item['key'],
                efnam=item['efnam'],
                ra=item['ra'],
                dec=item['dec'],
                epoch=item['epoch'],
                area_deg2=item['area_deg2'],
                fill_factor=item['fill_factor'],
                mag_lim=item['mag_lim'],
                obs_code=item['obs_code'],
                rate_mid_asphr=item['rate_mid_asphr'],
                _lib=lib,
            )
            by_index.append(p)
            pointings[pid] = p

        return cls(
            directory,
            pointings,
            by_index,
            n_epochs=int(n_epochs),
            detections_required=int(detections_required),
        )

    def keys(self) -> List[str]:
        return list(self.pointings.keys())

    def __getitem__(self, key: Union[str, int]) -> Pointing:
        if isinstance(key, int):
            return self.by_index[key]
        return self.pointings[key]

    def __iter__(self):
        return iter(self.by_index)

    def __len__(self) -> int:
        return len(self.by_index)

    def coverage_vs_magnitude(
        self,
        mags: MagLike,
        rate_asphr: Optional[float] = None,
        pointing_ids: Optional[Iterable[str]] = None,
        include_fill: bool = True,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Sum effective area over pointings vs magnitude.

        Args:
            mags: magnitude grid.
            rate_asphr: on-sky rate ["/hr]; default is each pointing's rate_cut midpoint.
            pointing_ids: Subset of pointing keys; default is all pointings.
            include_fill: multiply by fill factor when True.

        Returns:
            (mags array, total effective area array) in square degrees.
        """
        m = np.asarray(mags, dtype=float)
        if pointing_ids is None:
            selected = list(self.by_index)
        else:
            selected = [self.pointings[k] for k in pointing_ids]
        total = np.zeros(m.shape, dtype=float)
        for p in selected:
            total = total + np.asarray(
                p.effective_area(m, rate_asphr=rate_asphr, include_fill=include_fill),
                dtype=float,
            )
        return m, total
