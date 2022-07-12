"""
The OSSOS Survey Simulator module.  This is front-end for the OSSSim.  The primary OSSSim is writen in f95.
"""
import random

import numpy as np
from astropy import units as u
from astropy.time import Time
from astropy.units import Quantity
import SurveySubsF95

import ssim_util

T_ORB_M_UNITS = {'a': u.au,
                 'e': u.dimensionless_unscaled,
                 'inc': u.radian,
                 'node': u.radian,
                 'peri': u.radian,
                 'm': u.radian,
                 'epoch': u.day
                 }


class OSSSIM(object):
    """
    Outer Solar System Survey Simulator.

    This class simulates the process of observing a model of the solar system using a set of characterize observations.

    The model is provided to the class as an iterator.
    The survey characterization is provided as a SSimConf object with holds location and info about of characterization files

    :model:
        model iterator must provide a dictionary like object that provides Quantities for the following properties.
            - 'a', 'e', 'inc', 'node', 'peri', 'm', 'h'

        The iterator can also provide, optionally, Quantities for LIGHT_CURVE effects:
            - 'epoch', 'gb', 'period', 'amp', 'phase'

        These additional values will override the default values from light_curve_param values set in the __init__

        The iterator can also provide a color dictionary for types of objects being simulated, each object can have its own color
        dictionary or a global set can be provided.
            - 'colors'

        e.g. COLORS = {'g-x': 0.0,
              'r-x': -0.70,
              'i-x': -1.2,
              'z-x': -1.7,
              'u-x': 0.8,
              'V-x': 0.5,
              'B-x': 0.1,
              'R-x': -0.8,
              'I-x': -1.2,
              'x-x': 0.0}

    :characterization:
        See SSimConf documentation for details.

    """

    # default lighcurve parameters, appropriate for typical TNO
    LIGHT_CURVE_PARAMS = {
        'phase': 0.0,
        'period': 0.0 * u.hour,
        'gb': -0.12,
        'amplitude': 0.0 * u.mag,
    }

    def __init__(self, model, characterization_directory,
                 light_curve_params=None):
        """

        :param model: iterator that provides a dictionary of orbital elements and value of 'h'
        :type model: ssim_util.SSimModelFile
        :param characterization_directory: the path to survey characterization to be used.
        :type characterization_directory: str

        The colors dictionary declares the color of the KBO in multiple filters.  The 'x' refers to the bandpass that the values of 'H'
        (provided by the model iterator) are provided in.  Thus, if the model provide H_r value the colors given should be such that
        'x' is replaced with 'r'.  The default vector expects that the "H" values are given in 'g' filter and colors are fixed.

        The light_curve_params dictionary is a set of parameters to declare the phase-effect and light-curve variation in the flux.

        """
        self._colors = None
        if light_curve_params is None:
            light_curve_params = OSSSIM.LIGHT_CURVE_PARAMS
        self.light_curve_params = light_curve_params

        self.model = model
        self.characterization_directory = characterization_directory

    def simulate(self, detect_filename, seed=None, ntrack=None):
        """
        :param detect_filename: name of file to store dectections
        :type detect_filename: str
        :param seed: random seed passed to fortran modules, to allow reproducible simulations.
        :type seed: int
        :param ntrack: number of objects to track:
        :type ntrack: int

        ntrack has special values: 0 - till model exhausted, -ve till this many iterations, +ve till this many detections

        The detections will be listed in detect_filename which will have the columns defined in ssim_util.DetectFile.colnames
            units: degree, au, jd  (RA in hours, DEC in degrees)

        """

        # if no seed is provide generate a 'random' one... this is the seed passed to the fortran code.
        if seed is None:
            seed = random.randint(0, 123456789)

        detect_file = ssim_util.DetectFile(detect_filename)
        detect_file.epoch = self.model.epoch.jd
        detect_file.lambda_neptune = self.model.lambda_neptune
        detect_file.colors = self.model.colors
        detect_file.write_header(seed)

        n_iter = n_hits = n_track = 0

        for row in self.model:
            # pack the orbit into a t_orb_m object to pass to fortran module.
            o_m = SurveySubsF95.datadec.t_orb_m()
            for colname in row:
                if hasattr(o_m, colname):
                    if isinstance(row[colname], Quantity):
                        value = row[colname].to(T_ORB_M_UNITS[colname]).value
                    else:
                        value = row[colname]
                    setattr(o_m, colname, value)

            # attempt to detect this object
            # colors = self.colors
            colors = row.get('colors', self.model.colors)
            epoch = row.get('epoch', self.model.epoch.jd)
            row['gb'] = row.get('gb', self.light_curve_params['gb'])
            row['phase'] = row.get('phase', self.light_curve_params['phase'])
            row['period'] = row.get('period', self.light_curve_params['period'])
            row['amplitude'] = row.get('amplitude', self.light_curve_params['amplitude'])

            row['flag'], row['RA'], row['DEC'], row['d_ra'], row['d_dec'], row['r'], row['delta'], \
                row['m_int'], row['m_rand'], row['eff'], isur, row['M'], jdayp, ic, row['Survey'], \
                row['h_rand'], ierr = \
                SurveySubsF95.Surveysub.detos1(o_m,
                                               epoch,
                                               row['H'].to(u.mag).value,
                                               colors,
                                               row['gb'],
                                               row['phase'],
                                               row['period'].to(u.hour).value,
                                               row['amplitude'].to(u.mag).value,
                                               self.characterization_directory,
                                               seed)
            if ierr != 0:
                raise IOError("SSim failed with error code: {ierr}")

            # record model objects that are detected include flag if 'tracked'
            if row['flag'] > 0:
                # ic gives the filter that the object was 'detected' in, this allows us to determine the color
                row['color'] = colors[ic - 1]
                row['q'] = row['a'] * (1 - row['e'])
                row['M'] *= u.rad
                row['RA'] *= u.rad
                row['DEC'] *= u.rad
                row['Survey'] = row['Survey'].decode('utf-8')
                # m_int and h are in "x" band (filter of object creation)
                # m_rand and h_rand are in discovery filter
                n_hits += 1
                detect_file.write_row(row)
                if (row['flag'] > 2) and (np.mod(row['flag'], 2) == 0):
                    n_track += 1

            # stop the loop when maximum/desired detections have occurred if ntrack==0 then go to end of model file.
            if (0 < ntrack <= n_track) | (0 < -ntrack <= n_iter):
                break
        detect_file.write_footer(n_iter=n_iter, n_hits=n_hits, n_track=n_track)
