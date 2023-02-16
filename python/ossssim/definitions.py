"""
Define various default parameters and units.
"""
from collections import OrderedDict
from astropy import units

# default lightcurve parameters, appropriate for typical TNO
LIGHT_CURVE_PARAMS = {
    'phase': 0.0 * units.rad,
    'period': 0.0 * units.hour,
    'gb': -0.12 * units.mag,
    'amplitude': 0.0 * units.mag,
}

u = units


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



# values here are from the CFEPS L4j block
# 0.025 is the minimum unceratainty
# 0.272 is the uncertainty at 24.3
# 0.2 is the slope of mag getting brigter faintward of 23.8
# -0.2 is the slope of mag getting fainter  faintward of 23.8
MAG_ERR_PARAMS = [0.025, 0.272, 24.3, 0.2, 23.8,-0.2 ]

# assumes x == g
COLORS = OrderedDict(
    (('g-x', +0.0 * units.mag),
     ('r-x', -0.7 * units.mag),
     ('i-x', -1.2 * units.mag),
     ('z-x', -1.7 * units.mag),
     ('u-x', +0.3 * units.mag),
     ('V-x', +0.5 * units.mag),
     ('B-x', +0.1 * units.mag),
     ('R-x', -0.8 * units.mag),
     ('I-x', -1.2 * units.mag),
     ('x-x', +0.0 * units.mag)))

COLUMN_MAP = {'i': 'inc'}

colunits = {'a': units.au,
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
            }

COLUMN_WIDTH = 10

formats = {'flag': f'{COLUMN_WIDTH}d',
           'Survey': f'{COLUMN_WIDTH}s',
           'Comments': f'{COLUMN_WIDTH}s',
           'comp': f'{COLUMN_WIDTH}s',
           'RA': f'{COLUMN_WIDTH}.5f',
           'DEC': f'{COLUMN_WIDTH}.4f',
           'j': f'{COLUMN_WIDTH}',
           'k': f'{COLUMN_WIDTH}',
           'default': f'{COLUMN_WIDTH}.4f'}

Neptune = {'Longitude': 333.5078 * units.deg,
           'Epoch': 2456293.5 * units.day
           }
