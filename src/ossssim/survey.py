"""
Survey characterization access: areas, fill factors, and efficiency vs depth.

Example
-------
>>> from ossssim import SurveyCharacterization, Characterizations
>>> import numpy as np
>>> survey = SurveyCharacterization.from_directory(
...     Characterizations.surveys['CFEPS'])
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
    """One survey pointing with geometric area, fill factor, and efficiency."""

    id: str
    index: int  # 1-based Fortran index into the loaded survey
    efnam: str
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
            eta, _maglim = lib.pointing_eta(self.index, float(arr), rate)
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


class SurveyCharacterization:
    """
    Loaded survey characterization directory with keyed pointing access.

    Pointing keys use the efficiency-file basename without ``.eff``. When
    several pointings share the same efficiency file they are disambiguated
    as ``name#0``, ``name#1``, ... in ``pointings.list`` order. Use
    ``by_index`` for unambiguous iteration.
    """

    def __init__(
        self,
        directory: Path,
        pointings: Dict[str, Pointing],
        by_index: List[Pointing],
    ):
        self.directory = Path(directory)
        self.pointings = pointings
        self.by_index = by_index

    @classmethod
    def from_directory(
        cls,
        directory: Union[str, Path],
        lun: int = 21,
    ) -> "SurveyCharacterization":
        """Load ``pointings.list`` and associated ``.eff`` files via Fortran."""
        directory = Path(directory).resolve()
        if not (directory / 'pointings.list').is_file():
            raise FileNotFoundError(f"No pointings.list in {directory}")

        lib = ossssimlib.surveysub
        lib.reset_simulator()
        n_sur, ierr = lib.survey_load(str(directory), lun)
        if ierr != 0 or n_sur <= 0:
            raise RuntimeError(
                f"survey_load failed for {directory}: ierr={ierr}, n_sur={n_sur}"
            )

        # First pass: collect raw metadata and count duplicate efnam stems
        raw = []
        counts: Dict[str, int] = {}
        for i in range(1, n_sur + 1):
            area, fill = lib.pointing_geom(i)
            ra, dec = lib.pointing_center(i)
            efnam, epoch, code, mag_lim, rate_mid = lib.pointing_meta(i)
            efnam_s = _as_str(efnam)
            stem = Path(efnam_s).name
            if stem.lower().endswith('.eff'):
                stem = stem[:-4]
            counts[stem] = counts.get(stem, 0) + 1
            raw.append(
                dict(
                    index=i,
                    efnam=efnam_s,
                    stem=stem,
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

        # Second pass: assign ids (plain stem if unique, else stem#k)
        seen: Dict[str, int] = {}
        by_index: List[Pointing] = []
        pointings: Dict[str, Pointing] = {}
        for item in raw:
            stem = item['stem']
            if counts[stem] == 1:
                pid = stem
            else:
                k = seen.get(stem, 0)
                seen[stem] = k + 1
                pid = f"{stem}#{k}"
            p = Pointing(
                id=pid,
                index=item['index'],
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

        return cls(directory, pointings, by_index)

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
        Effective survey area as a function of magnitude.

        A_eff(m) = sum_i area_i * fill_i * eta_i(m, rate_i)

        Parameters
        ----------
        mags :
            Magnitude grid.
        rate_asphr :
            On-sky rate ["/hr] applied to every pointing. When None, each
            pointing uses its own rate_cut midpoint.
        pointing_ids :
            Subset of pointing keys; default is all pointings.
        include_fill :
            If False, omit fill factors from the sum.

        Returns
        -------
        mags, A_eff : ndarray
            Magnitude grid and effective area [deg^2].
        """
        m = np.asarray(mags, dtype=float)
        if pointing_ids is None:
            selected = self.by_index
        else:
            selected = [self.pointings[k] for k in pointing_ids]

        total = np.zeros(m.shape, dtype=float)
        for p in selected:
            total = total + np.asarray(
                p.effective_area(m, rate_asphr=rate_asphr, include_fill=include_fill),
                dtype=float,
            )
        return m, total
