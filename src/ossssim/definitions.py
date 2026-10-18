"""
Define various default parameters and units.
"""
from collections import OrderedDict

import numpy
from astropy import units
from astropy.time import Time

# default lightcurve parameters, appropriate for typical TNO
LIGHT_CURVE_PARAMS = {
    'phase': 0.0 * units.rad,
    'period': 0.0 * units.hour,
    'gb': -0.12 * units.mag,
    'amplitude': 0.0 * units.mag,
}

u = units

DEFAULT_SEED = 123456789
H_MIN = 1
H_MAX = 9
DEFAULT_SIZE = 1000
DEFAULT_MODEL_BAND = 'r'

# these are the units that the SurveySubF95 module (ie the fortran code) expects elements to be in.
T_ORB_M_UNITS = {'a': u.au,
                 'e': u.dimensionless_unscaled,
                 'inc': u.radian,
                 'node': u.radian,
                 'Node': u.radian,
                 'peri': u.radian,
                 'M': u.radian,
                 'epoch': u.day,
                 'x': u.au,
                 'y': u.au,
                 'z': u.au,
                 'vx': u.au/u.year,
                 'vy': u.au/u.year,
                 'vz': u.au/u.year
                 }



# values here are from the CFEPS L4j block
# 0.025 is the minimum uncertainty
# 0.272 is the uncertainty at 24.3
# 0.2 is the slope of mag getting brighter faintward of 23.8
# -0.2 is the slope of mag getting fainter  faintward of 23.8
MAG_ERR_PARAMS = [0.025, 0.272, 24.3, 0.2, 23.8, -0.2]


COLUMN_MAP = {'i': 'inc'}


column_description = {
    'a': 'semi-major axis',
    'e': 'eccentricity',
    'inc': 'inclination',
    'q': 'perihelion distance',
    'r': 'distance from sun',
    'M': 'mean anomaly at epoch of model',
    'node': 'longitude of ascending node',
    'peri': 'argument of perihelion',
    'Mt': 'Mean anomaly at Epoch of detection',
    'm_rand': 'magnitude in survey filter band sampled from uncertainty range',
    'h_rand': 'H magnitude in survey filter band sampled from uncertainty range',
    'color': 'model color of model object',
    'flag': '0==not detected, 1==detected, 2==tracked',
    'delta': 'distance from observer at detection',
    'm_int': 'intrinsic magnitude in filter band of the model',
    'h': 'absolute magnitude in filter band of the model',
    'H': 'absolute magnitude in filter band of the model',
    'eff': 'detection efficiency of field detected in',
    'RA': 'right ascension at detection',
    'DEC': 'declination at detection',
    'Survey': 'detection key survey/block (/ is a key delimiter, not a path)',
    'comp': 'component of orbital population model object is from',
    'dist': 'distance from observer',
    'j': 'j component of Neptune j:k exterior resonance',
    'k': 'k component of Neptune j:k exterior resonance',
    'phi': 'phase angle at detection',
    'resamp': 'amplitude of resonance libration',
    'Name': 'name of object',
    'n': 'mean motion',
    'Q': 'aphelion distance',
    'P': 'orbital period',
    'epoch': 'time at detection',
    'Epoch': 'time at detection',
    'x': 'Heliocentric x position at epoch of detection',
    'y': 'Heliocentric y position at epoch of detection',
    'z': 'Heliocentric z position at epoch of detection',
    'vx': 'Heliocentric x velocity at epoch of detection',
    'vy': 'Heliocentric y velocity at epoch of detection',
    'vz': 'Heliocentric z velocity at epoch of detection',
    'delta_v': 'velocity difference between object and observer',
    'dt': 'time between detection and epoch of orbit',
    'band': 'Bandpass of survey the object was detected in'
}

column_unit = {
    'a': units.au,
    'e': None,
    'inc': units.degree,
    'q': units.au,
    'r': units.au,
    'M': units.degree,
    'node': units.degree,
    'peri': units.degree,
    'Mt': units.degree,
    'm_rand': units.mag,
    'h_rand': units.mag,
    'color': units.mag,
     'flag': None,
    'delta': units.au,
    'm_int': units.mag,
    'h': units.mag,
    'H': units.mag,
    'eff': None,
    'RA': units.deg,
    'DEC': units.deg,
    'Survey': None,
     'comp': None,
    'dist': units.au,
    'j': None,
    'k': None,
    'phi': units.deg,
    'resamp': units.deg,
    'Name': None,
    'n': None,
    'Q': units.au,
    'P': units.year,
    'epoch': units.day,
    'Epoch': units.day,
    'x': u.au,
    'y': u.au,
    'z': u.au,
    'vx': u.au / u.year,
    'vy': u.au / u.year,
    'vz': u.au / u.year,
    'delta_v': u.meter / u.second,
    'dt': u.year,
    'band': None,
}

observables = ['RA', 'DEC', 'd_ra', 'd_dec', 'r', 'delta', 'm_int', 'Survey', 'eff',
               'm_rand', 'h_rand', 'Mt', 'h_rand', 'x', 'y', 'z', 'vx', 'vy', 'vz', 'band', 'color']

COLUMN_WIDTH = 14
# Detection key "survey/block" (name_len + 1 + block_len in Fortran parameters).
SURVEY_KEY_WIDTH = 65

column_format = {
    'flag': f'{COLUMN_WIDTH}d',
    'Survey': f'{SURVEY_KEY_WIDTH}s',
    'Comments': f'{COLUMN_WIDTH}s',
    'comp': f'{COLUMN_WIDTH}s',
    'RA': f'>{COLUMN_WIDTH}.5f',
    'DEC': f'>{COLUMN_WIDTH}.4f',
    'j': f'{COLUMN_WIDTH}',
    'k': f'{COLUMN_WIDTH}',
    'default': f' >{COLUMN_WIDTH}.5f',
    'band': f'{COLUMN_WIDTH}s'
}
for col in column_description:
    column_format[col] = column_format.get(col, column_format['default'])

column_dtype = {
    'flag': 'i',
    'Survey': 'U',
    'Comments': 'U',
    'comp': 'U',
    'RA': 'f',
    'DEC': 'f',
    'j': 'i',
    'k': 'i',
    'default': 'f',
    'band': 'U'
}

Neptune = {
    'Longitude': 332.3354 * units.deg,
    'Latitude': -0.6214 * units.deg,
    'Distance': 29.99064931353 * units.au,
    'Speed': -0.0552043 * units.km / units.s,
    'Epoch': Time(2456293.5, format='jd'),
    'SemimajorAxis': 30.06952752 * units.au,
}
