"""
The OSSOS Survey Simulator module.  This is front-end for the OSSSim.  The primary OSSSim is writen in f95.
The OSSOS Survey Simulator module.  This is front-end for the OSSSim.  The primary OSSSim is writen in f95.
"""
import random
from astropy import units as u
from astropy.units import Quantity
from astropy.table import Row
from .lib import SurveySubsF95
from . import definitions

# these are the units that the SurveySubF95 module (ie the fortran code) expects elements to be in.
T_ORB_M_UNITS = {'a': u.au,
                 'e': u.dimensionless_unscaled,
                 'inc': u.radian,
                 'node': u.radian,
                 'Node': u.radian,
                 'peri': u.radian,
                 'M': u.radian,
                 'epoch': u.day
                 }


class OSSSSim:
    """
    Outer Solar System Survey Simulator.

    This class simulates the process of observing a model of the solar system using a set of characterized observations.
    """


    def __init__(self, characterization_directory):
        """
        Args:
            characterization_directory (str): the path to survey characterization to be used.

        Format of the characterization_directory is described at https://github.com/OSSOS/SurveySimulator/tree/master/Surveys

        """

        self.characterization_directory = characterization_directory

    def simulate(self, row, seed=None, epoch=None, colors=None):
        """
        Pass the target elements to detos1 and determine if this target would be detected.

        Args:
            row (Row or dict): elements of the target to simulate
            seed (int): a seed to pass to the fortran code to allow reproducible simulations
            epoch (Quantity or float): JD of elements, can also be provided for each row sent to simulate
            colors (dict): colors of the target index by 'g-x' etc, see below.

        Returns:
            Row or dict: the target elements/values at time of simulated detection.

        row should have values for a, e, inc, node, peri, M, H, [epoch]
        can, optionally, also define: gb, phase, period, amplitude

        the result row has, in addition to above, the following items:

        flag: 0 - not detected, 1 - detected, 2 - tracked

        following are values at time of detection, None if not detected

        RA: RA of target

        DEC: DEC of target

        d_ra: rate of RA sky motion

        d_dec: rate of DEC sky motion

        delta: distance from Earth

        r: distance from Sun

        m_int: intrinsic magnitude in filter of target H (ie. the model filter)

        m_rand: magnitude at detection, includes scatter due to flux measurement uncertainty (ie. the survey filter)

        h_rand: inferred absolute magnitude based on m_rand and detection circumstances

        eff: the detection field's efficiency of detection for a source of m_rand

        M: the mean anomaly at detection

        Survey: a string indicating which field detected the target

        The colors list declares the color of the KBO in multiple filters.  The 'x' refers to the bandpass that the values of 'H'
        (provided to simulate) are provided in.  Thus, if the target passed to simulate provides H_r value the colors given should
        be such that 'x' is replaced with 'r'.  The default vector expects that the "H" values are given in 'g' filter and colors
        are fixed.  See ModelFile for more details.

        A default light_curve_params dictionary is used if target provided to simulate doesn't have any, set in
        ossssim.definitions.LIGHT_CURVE_PARAMS
        """

        # if no seed is provide generate a 'random' one... this is the seed passed to the fortran code.
        if seed is None:
            seed = random.randint(0, 123456789)

        # pack the orbit into a t_orb_m object to pass to fortran module.
        o_m = SurveySubsF95.datadec.t_orb_m()
        row = dict(row)
        for colname in row:
            if hasattr(o_m, colname.lower()):
                if isinstance(row[colname], Quantity):
                    value = row[colname].to(T_ORB_M_UNITS[colname]).value
                else:
                    value = row[colname]
                setattr(o_m, colname.lower(), value)

        # attempt to detect this object
        if 'colors' in row.keys():
            colors = [x.to(u.mag).value for x in row['colors']]
        else:
            colors = [definitions.COLORS[x].to(u.mag).value for x in definitions.COLORS]
        # if isinstance(colors, dict):
        #    # order the list correctly using the keys returned by the OrderedDict ModelFile.COLORS
        #    colors = [ colors[x] for x in ModelFile.COLORS]
        # colors = [ x.to(u.mag).value for x in colors]

        epoch = 'epoch' in row.keys() and row['epoch'].to(u.day).value or epoch.to(u.day).value
        gb = 'gb' in row.keys() and row['gb'].to(u.mag).value or definitions.LIGHT_CURVE_PARAMS['gb'].to(u.mag).value
        phase = 'phase' in row.keys() and row['phase'].to(u.radian).value or definitions.LIGHT_CURVE_PARAMS['phase'].to(u.radian).value
        period = 'period' in row.keys() and row['period'].to(u.day).value or definitions.LIGHT_CURVE_PARAMS['period'].to(u.day).value
        amplitude = 'amplitude' in row.keys() and row.get('amplitude').to(u.mag).value or definitions.LIGHT_CURVE_PARAMS['amplitude'].to(u.mag).value
        H = 'H' in row.keys() and row['H'].to(u.mag).value or ('h' in row.keys() and row['h'].to(u.mag).value or 0.0)

        row['flag'], row['RA'], row['DEC'], row['d_ra'], row['d_dec'], row['r'], row['delta'], \
            row['m_int'], row['m_rand'], row['eff'], isur, row['Mt'], jdayp, ic, row['Survey'], \
            row['h_rand'], ierr = \
            SurveySubsF95.Surveysub.detos1(o_m,
                                           epoch,
                                           H,
                                           colors,
                                           gb,
                                           phase,
                                           period,
                                           amplitude,
                                           self.characterization_directory,
                                           seed)
        if ierr != 0:
            raise IOError(f"SSim failed with error code: {ierr}")

        row['delta'] *= u.au
        row['r'] *= u.au
        row['m_int'] *= u.mag
        row['m_rand'] *= u.mag
        row['h_rand'] *= u.mag
        # ic gives the filter that the object was 'detected' in,
        # this allows us to determine the color of target
        row['color'] = colors[ic - 1] * u.mag
        row['q'] = row['a'] * (1 - row['e'])
        row['Mt'] *= u.rad
        row['RA'] *= u.rad
        row['DEC'] *= u.rad
        row['Survey'] = row['Survey'].decode('utf-8')

        # m_int and h are in "x" band (filter of object creation)
        # m_rand and h_rand are in discovery filter
        return row
