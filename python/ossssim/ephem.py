import math

from astropy.coordinates import SkyCoord
from astropy.time import Time
from astropy.units import Quantity

from .lib import SurveySubsF95
from . import definitions

T_ORB_M_UNITS = definitions.T_ORB_M_UNITS


class Ephem:
    """
    Provides the coordinate of an object given the observing code, julian date and elements.
    """

    def __init__(self, epoch:float, code:int = 500):
        """
        param row: a dictionary holding a row of elements from a ModelFile object
        param code: observatory code to compute observed location relative to.
        """
        self.epoch = epoch
        self.code = code
        self._o_m = None
        self._pos = None
        self._obspos = None
        self._coord = None
        self.row = None

    def __call__(self, row):
        self.row = row
        self._o_m = None
        self._coord = None
        return self.coord

    @property
    def o_m(self):
        """
        Take a row of elements and pack into form suitable to be sent to the F95 routines.
        """
        if self._o_m is None:
            self._o_m = SurveySubsF95.datadec.t_orb_m()
            row = dict(self.row)
            for colname in row:
                if hasattr(self._o_m, colname.lower()):
                    if isinstance(row[colname], Quantity):
                        value = row[colname].to(T_ORB_M_UNITS[colname]).value
                    else:
                        value = row[colname]
                    setattr(self._o_m, colname.lower(), value)
        return self._o_m

    @property
    def obspos(self) -> list:
        if self._obspos is None:
            self._obspos, vel, self._r, ierr = SurveySubsF95.xvutils.obspos(self.code,
                                                                            self.epoch)
        return self._obspos

    @property
    def coord(self) -> SkyCoord:
        if self._coord is None:
            pos = SurveySubsF95.elemutils.pos_cart(self.o_m)
            r = math.sqrt(pos.x**2 + pos.y**2 + pos.z**2)
            delta, ra, dec = SurveySubsF95.numutils.radececlxv(pos, self.obspos)
            return SkyCoord(ra,
                            dec,
                            distance=r,
                            unit=('radian', 'radian', 'au'),
                            obstime=Time(self.epoch, format='jd'))
