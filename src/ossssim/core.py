"""
The OSSOS Survey Simulator module.  This is front-end for the OSSSim.  The primary OSSSim is writen in f95.
The OSSOS Survey Simulator module.  This is front-end for the OSSSim.  The primary OSSSim is writen in f95.
"""
import logging
from astropy import units as u
from astropy.units import Quantity
from astropy.table import Row
from astropy.time import Time

from .color import PhotSpec
import ossssimlib
import rebound
import os
import numpy
from . import definitions
T_ORB_M_UNITS = definitions.T_ORB_M_UNITS


class Cartesian(object):
    """
    A class to provide the cartesian state vector of a set of objects.
    """
    def __init__(self, epoch=definitions.Neptune['Epoch']):
        self.epoch = epoch
        self.rebound_archive = f"Rebound_Archive{epoch.jd}.bin"
        self.initialize_sim()

    def _make_rebound_archive_file(self):
        _sim = rebound.Simulation()
        _sim.add("Sun",  date=self.epoch.to_datetime())
        _sim.add("Jupiter", date=self.epoch.to_datetime())
        _sim.add("Saturn", date=self.epoch.to_datetime())
        _sim.add("Uranus", date=self.epoch.to_datetime())
        _sim.add("Neptune", date=self.epoch.to_datetime())
        _sim.move_to_com()
        _sim.save_to_file(self.rebound_archive)
        del _sim

    def initialize_sim(self):
        """
        Start/open a rebound Simulation object used to compute the cartesian locations of the particles in this model.
        """
        if not os.access(self.rebound_archive, os.F_OK):
            self._make_rebound_archive_file()
        return rebound.Simulation(self.rebound_archive)

    def __call__(self, **kwargs) -> dict:
        """
        Provide the state vector of the orbits at Mt or M if Mt is not provided.


        Returns:
            dictionary of x/y/z/vx/vy/vz
        """
        sim = self.initialize_sim()
        if "rows" not in kwargs:
            return self.single_row(**kwargs)

        rows = kwargs['rows']
        M = 'Mt' in rows and 'Mt' or 'M'
        for idx in range(len(rows['a'])):
            sim.add(a=rows['a'][idx].to('au').value,
                    e=rows['e'][idx],
                    inc=rows['inc'][idx].to('rad').value,
                    Omega=rows['node'][idx].to('rad').value,
                    omega=rows['peri'][idx].to('rad').value,
                    M=rows[M][idx].to('rad').value,
                    )

        for s in ['x', 'y', 'z', 'vx', 'vy', 'vz']:
            rows[s] = [p.__getattribute__(s) * definitions.column_unit[s] for p in sim.particles[5:]]
        return rows

    def single_row(self, **kwargs):
        """
        Provide the state vector of the orbits at Mt or M if Mt is not provided.

        Returns:
            dictionary of x/y/z/vx/vy/vz
        """
        sim = self.initialize_sim()
        sim.add(a=kwargs['a'].to('au').value,
                e=kwargs['e'],
                inc=kwargs['inc'].to('rad').value,
                Omega=kwargs['node'].to('rad').value,
                omega=kwargs['peri'].to('rad').value,
                M=kwargs.get('Mt', kwargs.get('M')).to('rad').value
                )
        cartesian = {}
        for s in ['x', 'y', 'z']:
            cartesian[s] = sim.particles[-1].__getattribute__(s) * u.au
        for s in ['vx', 'vy', 'vz']:
            cartesian[s] = sim.particles[-1].__getattribute__(s) * (2 * numpy.pi) * u.au / u.year
        return cartesian


class OSSSSim:
    """
    Outer Solar System Survey Simulator.

    This class simulates the process of observing a model of the solar system using a set of characterized observations.
    """

    def __init__(self, characterization_directory, seed):
        """
        Args:
            characterization_directory (str): path to a survey characterization directory.
                Full surveys are distributed separately in SurveySimulator-Data;
                format docs are under docs/; small fixtures live in F95/tests/Surveys/.

        """

        self.characterization_directory = characterization_directory
        self.cartesian = Cartesian(epoch=definitions.Neptune['Epoch'])
        self.seed = seed 
        ossssimlib.surveysub.reset_simulator()


    def simulate(self, row: dict, colors: PhotSpec, model_band, seed=None, epoch=None, debug=False) -> dict:
        """
        Pass the target elements to detos1 and determine if this target would be detected.

        Args:
            row (Row or dict): elements of the target to simulate
            colors (PhotSpec): colors of the target index by 'g-x' etc, see below.
            model_band (str): a single letter designator of the bandpass that the H values are provided in.
            seed (int): a seed to pass to the fortran code to allow reproducible simulations
            epoch (Quantity or float or Time): JD of elements, can also be provided for each row sent to simulate
            debug (bool): if True, have Fortran detos routine create a log of detection process.

        Returns:
            Row or dict: the target elements/values at time of simulated detection.

        row should have values for a, e, inc, node, peri, M, H, [epoch]
        can, optionally, also define: gb, phase, period, amplitude

        the result row has, in addition to above, the following items:

        flag: 0 - not detected, 1 - detected, 2 - tracked, 3 - characterized and lost, 4 - tacked and characterised

        When comparing to a detected sample you likely want to only use those with flag == 4

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

        The colors list declares the color of the KBO in multiple filters. The value of model_band_pass is
        used to transform the color dictionary into colors in the bandpass of the model.

        A default light_curve_params dictionary is used if target provided to simulate doesn't have any, set in
        ossssim.definitions.LIGHT_CURVE_PARAMS
        """

        if seed is not None:
            logging.warning('seed is now set in the constructor. ignored here')

        # pack the orbit into a t_orb_m object to pass to fortran module.
        o_m = ossssimlib.datadec.t_orb_m()
        row = dict(row)
        for colname in row:
            if hasattr(o_m, colname.lower()):
                if isinstance(row[colname], Quantity):
                    value = row[colname].to(T_ORB_M_UNITS[colname]).value
                else:
                    value = row[colname]
                setattr(o_m, colname.lower(), value)

        # attempt to detect this object
        color_offset_array = colors.colors_list(row['comp'], model_band)
        if 'epoch' in row.keys():
            epoch_jd = row['epoch']
        else:
            epoch_jd = epoch
        if isinstance(epoch_jd, Quantity):
            epoch_jd = epoch_jd.to(u.day).value
        if isinstance(epoch_jd, Time):
            epoch_jd = epoch_jd.jd
        gb = 'gb' in row.keys() and row['gb'].to(u.mag).value or definitions.LIGHT_CURVE_PARAMS['gb'].to(u.mag).value
        phase = 'phase' in row.keys() and row['phase'].to(u.radian).value or definitions.LIGHT_CURVE_PARAMS['phase'].to(u.radian).value
        period = 'period' in row.keys() and row['period'].to(u.day).value or definitions.LIGHT_CURVE_PARAMS['period'].to(u.day).value
        amplitude = 'amplitude' in row.keys() and row.get('amplitude').to(u.mag).value or definitions.LIGHT_CURVE_PARAMS['amplitude'].to(u.mag).value
        H = 'H' in row.keys() and row['H'].to(u.mag).value or ('h' in row.keys() and row['h'].to(u.mag).value or None)
        if H is None:
            raise ValueError("H or h must be provided to simulate")

        row['flag'], row['RA'], row['DEC'], row['d_ra'], row['d_dec'], row['r'], row['delta'], \
            row['m_int'], row['m_rand'], row['eff'], isur, row['Mt'], jdayp, ic, row['Survey'], \
            row['h_rand'], ierr = \
            ossssimlib.surveysub.detos1(o_m,
                                           epoch_jd,
                                           H,
                                           color_offset_array,
                                           gb,
                                           phase,
                                           period,
                                           amplitude,
                                           self.characterization_directory,
                                           self.seed,
                                           debug)
        if ierr != 0:
            raise IOError(f"SSim failed with error code: {ierr}")
        row['j'] = row.get('j', 0)
        row['k'] = row.get('k', 0)
        row['delta'] *= u.au
        row['r'] *= u.au
        row['m_int'] *= u.mag
        row['m_rand'] *= u.mag
        row['h_rand'] *= u.mag
        # ic is the 1-based Fortran filter_to_index of the discovery bandpass.
        # Convert to band letter and look up the matching 0-based colors_list slot.
        if ic > 0:
            band = chr(ic + ord('A') - 1)
            row['band'] = band
            row['color'] = color_offset_array[ord(band) - ord('A')] * u.mag
        else:
            row['band'] = ''
            row['color'] = 0*u.mag
        row['q'] = row['a'] * (1 - row['e'])
        row['Mt'] *= u.rad
        row['RA'] *= u.rad
        row['DEC'] *= u.rad
        row['Survey'] = row['Survey'].decode('utf-8')
        if row['flag'] > 0:
            row.update(self.cartesian(**row))
        numpy.ma.array(row, mask=False)

        # m_int and h are in "x" band (filter of object creation)
        # m_rand and h_rand are in discovery filter
        return row
