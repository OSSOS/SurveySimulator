"""
Define various default parameters and units.
"""
from collections import OrderedDict
from astropy import units as u

COLORS = OrderedDict(
    (('g-x', +0.0 * u.mag),
     ('r-x', -0.7 * u.mag),
     ('i-x', -1.2 * u.mag),
     ('z-x', -1.7 * u.mag),
     ('u-x', +0.8 * u.mag),
     ('V-x', +0.5 * u.mag),
     ('B-x', +0.1 * u.mag),
     ('R-x', -0.8 * u.mag),
     ('I-x', -1.2 * u.mag),
     ('x-x', +0.0 * u.mag)))

COLUMN_MAP = {'i': 'inc'}

colunits = {'a': u.au,
            'e': None,
            'inc': u.degree,
            'q': u.au,
            'r': u.au,
            'm': u.degree,
            'node': u.degree,
            'peri': u.degree,
            'mt': u.degree,
            'm_rand': u.mag,
            'h_rand': u.mag,
            'color': u.mag,
            'flag': None,
            'delta': u.au,
            'm_int': u.mag,
            'h': u.mag,
            'H': u.mag,
            'eff': None,
            'RA': u.deg,
            'DEC': u.deg,
            'Survey': None,
            'comp': None,
            'dist': u.au,
            'j': None,
            'k': None
            }
