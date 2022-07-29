"""
Model classes for the Outer Solar System Survey Simulator.
"""
import logging
import os
import re
import time
from abc import ABC, abstractmethod
from collections import OrderedDict
from collections.abc import Iterable

import numpy
from astropy import units
from astropy.table import QTable
from astropy.time import Time
from astropy.units import Quantity
from numpy import random

from . import definitions
from . import distributions

RE_FLOAT = re.compile('\d+\.?\d*[de]?\d*')


def get_floats_in_str(line):
    """
    Scan a string (line) and return all float like strings as array of floats.

    Converts floats expressed as fortran doubles like '1.0d0' to '1.0e0' before conversion.
    """
    values = RE_FLOAT.findall(line)

    try:
        # convert values to floats if possible.
        result = [float(x.replace('d', 'e')) for x in values]
    except ValueError as ex:
        logging.debug(ex)
        logging.debug(f"Failed to parse: {line}")
        result = values
    return result


class ResultsFile:
    """
    ModelFile structure for output file from Simulator detections.
    """
    colnames = ['a', 'e', 'inc', 'node', 'peri', 'M', 'H',
                'q', 'r', 'm_rand', 'H_rand',
                'color', 'comp', 'j', 'k']

    def __init__(self, filename, randomize=False):
        self.filename = filename
        self.randomize = randomize
        self._header = None
        self._header_parsed = False
        self._colnames = None
        self._colors = None
        self._epoch = None
        self._longitude_neptune = None
        self._seed = None
        self.header_lines = []
        self.file_location = None

    @property
    def epoch(self):
        """
        Epoch of coordinates of orbit
        """
        if self._epoch is None:
            raise ValueError(f"epoch has not been set.")
        return self._epoch

    @epoch.setter
    def epoch(self, value):
        """
        Args:
            value (Time or float or int or Quantity): the value of the epoch
        """
        if isinstance(value, Time):
            self._epoch = value
        elif isinstance(value, float) or isinstance(value, int):
            self._epoch = Time(value, format='jd')
        elif isinstance(value, Quantity):
            self._epoch = Time(value.to('day').value, format='jd')
        else:
            raise ValueError(f"Don't know how to set epoch using: {value}")

    @property
    def longitude_neptune(self):
        """
        Longitude of Neptune at epoch
        """
        if self._longitude_neptune is None:
            raise ValueError(f"longitude_neptune has not been set.")
        return self._longitude_neptune

    @longitude_neptune.setter
    def longitude_neptune(self, value):
        if not isinstance(value, Quantity) and value is not None:
            raise ValueError(f"longitude_neptune must be set to a Quantity with units")
        self._longitude_neptune = value

    @property
    def colors(self):
        """Returns color array from file header or default if no color array in header."""
        if self._colors is None:
            # pickup the default values
            self._colors = list(definitions.COLORS.values())
        return self._colors

    @colors.setter
    def colors(self, values):
        self._colors = values

    def write_row(self, this_row):
        """
        Given a dictionary of row values write the row according to the order in colnames

        :param this_row: Dictionary of values to write to row.
        """
        with open(self.filename, 'a') as f_detect:
            # start sep as 2 spaces as column name header starts with '# '
            sep = "  "
            for colname in self.colnames:
                col_format = definitions.formats.get(colname, definitions.formats['default'])
                if isinstance(this_row[colname], Quantity):
                    value = this_row[colname].to(definitions.colunits[colname]).value
                else:
                    value = this_row[colname]
                o_str = "{sep}{value:{col_format}}".format(value=value,
                                                           col_format=col_format,
                                                           sep=sep)
                f_detect.write(o_str)
                sep = " "
            f_detect.write('\n')

    def write_header(self, seed):
        """
        :param seed: the seed for random generator that resulted in these detections.
        :type seed: int
        """
        with open(self.filename, 'w') as f_detect:
            f_detect.write(f"# Seed = {seed:10d}\n")
            f_detect.write(f"# Epoch of elements: JD = {self.epoch}\n")
            f_detect.write(f"# Longitude of Neptune: longitude_neptune = {self.longitude_neptune}\n")
            if self.colors is not None:
                color_str = [f"{c.to(units.mag).value:5.2f} " for c in self.colors]
                f_detect.write(f"# Colors = {color_str}\n")
            f_detect.write(f"#\n")
            date = time.strftime("%Y-%m-%dT%H:%M:%S.000  %z")
            f_detect.write(f"# Creation_time: {date:30s}\n")
            f_detect.write('# flag: >0: detected; >2: characterized; 0 mod(2): tracked\n')
            f_detect.write('# Survey: name of the observing block where model object was detected\n')
            f_detect.write("# ")
            for colname in self.colnames:
                f_detect.write("{colname:>{width}s} ".format(colname=colname,
                                                             width=definitions.COLUMN_WIDTH))
            f_detect.write("\n")

    def write_footer(self, n_iter, n_hits, n_track):
        """
        Write a footer with the results of the survey simulation

        This is done as a footer instead of header to allow streaming output.
        """
        with open(self.filename, 'a') as f_detect:
            f_detect.write('# Total number of objects:   {:11d}\n'.format(n_iter))
            f_detect.write('# Number of detections:      {:7d}\n'.format(n_hits))
            f_detect.write('# Number of tracked objects: {:7d}\n'.format(n_track))


class ModelOutputFile(ResultsFile):
    """
    Output format used to store the input model, used when model is parametric and you want to keep a record of input for diagnostics
    """

    colnames = ['a', 'e', 'inc', 'node', 'peri', 'M', 'H', 'q',
                'color', 'comp', 'j', 'k']


class DetectFile(ResultsFile):
    """
    Detected object output file structure.
    """
    colnames = ['a', 'e', 'inc', 'node', 'peri', 'M', 'H', 'q', 'r', 'Mt', 'm_rand', 'h_rand', 'color', 'flag',
                'delta', 'm_int', 'eff', 'RA', 'DEC', 'comp', 'j', 'k']


class TrackFile(ResultsFile):
    """
    Tracked object output file structure.
    """
    colnames = ['a', 'e', 'inc', 'node', 'peri', 'M', 'H', 'q', 'r',
                'Mt', 'm_rand', 'H_rand', 'color', 'Survey', 'comp', 'j', 'k']


class ModelFile(Iterable):
    """
    A class to drive the SSim using a standard model file.

    ModelFile opens file and reads the header for the epoch, seed, longitude_neptune and colors
    and then loops over or randomly offsets into the file to read model objects.
    """

    def __init__(self, filename, randomize=False):
        self.filename = filename
        self.randomize = randomize
        self._header = None
        self._header_parsed = False
        self._colnames = None
        self._colors = None
        self._epoch = None
        self._longitude_neptune = None
        self._seed = None
        self.header_lines = []
        self._f_obj = open(self.filename, 'r')
        self.f_loc = 0
        self._targets = None
        self._f = None

    @property
    def epoch(self):
        """
        Epoch of coordinates of orbit read from model file header.
        """
        if self._epoch is None:
            self._epoch = Time(float(self.header['JD'].replace('d', 'e')),
                               format='jd').jd * units.day
        return self._epoch

    @property
    def longitude_neptune(self):
        """
        Longitude of Neptune at epoch
        """
        if self._longitude_neptune is None:
            self._longitude_neptune = float(self.header['lambdaN'].replace('d', 'e')) * units.radian
        return self._longitude_neptune

    @property
    def colors(self):
        """Returns color array from file header or default if no color array in header."""
        if self._colors is None:
            # pickup the
            self._colors = list(definitions.COLORS.values())
            _header_colors = get_floats_in_str(self.header.get('colors', ""))
            for idx in range(len(_header_colors)):
                self._colors[idx] = _header_colors[idx]
        return self._colors

    @property
    def colnames(self):
        """
        Parse the file header (lines that start with #) and return the last header
        line split on spaces as the names of the data columns in filename.
        """
        if self._colnames is None:
            self._colnames = []
            for colname in self.header.get('colnames', "").split():
                colname = definitions.COLUMN_MAP.get(colname, colname)
                self._colnames.append(colname)
        if len(self._colnames) == 0:
            raise IOError(f"Failed to get column names in {self.filename}\n")
        return self._colnames

    @property
    def header(self):
        """
        Parse the first block of comment lines for header definition.
        """
        if self._header is not None or self._header_parsed:
            return self._header
        previous_line = None
        self._header = OrderedDict()
        with open(self.filename, 'r') as f_obj:
            if self.f_loc is not None:
                f_obj.seek(self.f_loc)
            self.f_loc = f_obj.tell()
            while True:
                line = f_obj.readline()
                # line = line.decode('utf-8')
                if line.startswith('#'):
                    logging.debug(f"Parsing Comment: {line}")
                    self.header_lines.append(line[1:])
                    if line.strip() == "#":
                        continue
                    if '=' in line:
                        keyword = line.split('=')[0].split()[-1]
                        value = line.split('=')[1].strip().split()[0]
                        self._header[keyword] = value
                    previous_line = line
                    self.f_loc = f_obj.tell()
                else:
                    if previous_line is not None:
                        # expect that the last header line is the column name header.
                        self._header['colnames'] = previous_line[1:]
                    break
            self._header_parsed = True
        if self._header is None:
            raise IOError(f"Failed to parse keywords from header of {self.filename}")
        self._header['_end_of_header_offset'] = self.f_loc
        return self._header

    def __iter__(self):
        return self

    def __next__(self):
        """
        Get the next line or a random line that is not a comment line from the file and parse into a row
        """
        if self.randomize:
            while True:
                # offset to random location in the file.
                self._f_obj.seek(random.randint(self.header['_end_of_header_offset'],
                                                os.stat(self.filename).st_size))
                try:
                    # read to the end of this line.
                    self._f_obj.readline()
                    while True:
                        line = self._f_obj.readline()
                        if len(line) == 0:
                            raise EOFError
                        if line[0] != "#":
                            break
                    break
                except EOFError:
                    self._f_obj.close()
                    self._f_obj = open(self.filename)
                    pass
        else:
            while True:
                line = self._f_obj.readline()
                if type(line) == bytes:
                    line = line.decode('utf-8')
                if len(line) == 0:
                    raise StopIteration
                if not line.startswith('#'):
                    break
        values = line.split()
        row = OrderedDict()
        for idx, colname in enumerate(self.colnames):
            try:
                if '.' in values[idx]:
                    value = float(values[idx].replace('d', 'e'))
                else:
                    value = int(values[idx])
            except ValueError:
                value = str(values[idx])
            except IndexError as ex:
                # for non-resonant we don't need to have j/k defined in file.
                if colname in ['j', 'k']:
                    value = 0
                else:
                    raise ex
            if definitions.colunits[colname] is not None:
                value = value * definitions.colunits[colname]
            row[colname] = value
        return row

    @property
    def targets(self):
        """
        targets set by looping over the entire file and returning a 'QTable'.  This can be used when you want access to the
        entire table of data rather than just reading one-line at a time.
        """

        if self._targets is None:
            # need to start from the top of the data range
            # get current file location.
            loc = self._f_obj.tell()

            # move the pointer to the end of the header.
            self._f_obj.seek(self.header['_end_of_header_offset'])
            values = {}
            for column in self.colnames:
                values[column] = []
            for row in self:
                for column in self.colnames:
                    values[column].append(row[column])
            self._targets = QTable(values)

            # set the file pointer back to where we were before.
            self._f_obj.seek(loc)

        return self._targets

    @property
    def f(self):
        """True anomaly of the orbit.  If the mean anomaly at detection (Mt) exists than use that
        else use the mean anomaly at model epoch

        Returns:
            (float): true anomaly
        """
        if self._f is None:
            e = self.targets['e']
            if 'M' not in self.targets.colnames:
                ValueError("No mean anomaly (M) in model?")
            if 'Mt' in self.targets.colnames:
                m = self.targets['Mt']
            else:
                m = self.targets['M']
            big_e = m - self.targets['e'] * numpy.sin(m) * units.rad
            converged = numpy.zeros(len(self.targets)) > 0
            f1 = numpy.zeros(len(self.targets)) * units.rad
            fp = numpy.ones(len(self.targets))
            n = 0
            while n < 1**6 and numpy.sum(~converged) > 0:
                f1[~converged] = big_e[~converged] - e[~converged] * numpy.sin(big_e[~converged]) * units.rad - m[~converged]
                f1[converged] = 0 * units.rad
                fp[~converged] = (1.0 - e[~converged] * numpy.sin(big_e[~converged]))
                fp[converged] = 1
                delta = -f1 / fp
                big_e = big_e + delta
                n += 1
                converged = delta < (1e-8 * units.rad)
            self._f = (2. * numpy.arctan(((1. + e) / (1. - e)) ** 0.5 * numpy.tan(big_e / 2.)))
        return self._f

    @property
    def cartesian(self):
        """
        provide the state vector of the orbits.

        Returns:
            (dict[ndarray, ndarray, ndarray, ndarray, ndarray, ndarray]): dictionary of x/y/z/vx/vy/vz
        """

        # compute the unit vector
        f = self.f
        capom = self.targets['node']
        om = self.targets['peri']
        inc = self.targets['inc']
        a = self.targets['a']
        e = self.targets['e']
        u = om + f
        np = numpy
        x_hat = np.cos(u) * np.cos(capom) - np.cos(inc) * np.sin(capom) * np.sin(u)
        y_hat = np.cos(u) * np.sin(capom) + np.cos(inc) * np.cos(capom) * np.sin(u)
        z_hat = np.sin(inc) * np.sin(u)

        # compute the angular momentum vector (unit vector)
        hx = np.sin(capom) * np.sin(inc)
        hy = -np.cos(capom) * np.sin(inc)
        hz = np.cos(inc)

        # assuming not parabolic, here the magnitudes of the vectors
        r = a * (1.0 - e * e) / (1.0 + e * np.cos(f))
        h = (a * (1.0 - e * e)) ** 0.5

        # position vectors
        x = r * x_hat
        y = r * y_hat
        z = r * z_hat

        # compute components of vector theta hat
        thx = hy * z_hat - hz * y_hat
        thy = hz * x_hat - hx * z_hat
        thz = hx * y_hat - hy * x_hat

        # obtain the velocity vector's components and calculate v
        theta_dot = h / (r * r)
        r_dot = e * np.sin(f) / h

        vx = r * theta_dot * thx + r_dot * x_hat
        vy = r * theta_dot * thy + r_dot * y_hat
        vz = r * theta_dot * thz + r_dot * z_hat

        return {'x': x, 'y': y, 'z': z, 'vx': vx, 'vy': vy, 'vz': vz}


class HDistribution:
    """
    Provide a class to describe the size distribution of the object being simulated.
    """
    def __init__(self, func, **kwargs):
        """
        Args:
            func (func): the function that will be used for the H-distribution.
            kwargs (dict): A dictionary of parameters used by the H distribution function.
        """
        self.params = kwargs
        self.func = func

    def __call__(self):
        return self.func(**self.params) * units.mag


class Parametric(ABC):
    """
    This abstract class defines methods needed to build a parametric Outer Solar System model for use as a model input for OSSSSim
    """

    def __init__(self,
                 size: int = 1000000,
                 seed: int = 123456789,
                 epoch: Quantity = 2456839.5 * units.day,
                 j: int = 0,
                 k: int = 0,
                 longitude_neptune: Quantity = 5.876 * units.rad,
                 **kwargs) -> None:
        """
        Setup the boundaries of the simulation.  size and seed are used to initialize a dist_utils.Distribution class.

        Args:
            size: Determines size of the arrays to be generated (default=10^6).
            seed: Initialize distributions with this seed to enable reproducible models.
            epoch: JD epoch of elements and Neptune longitude
            j: MMR Neptune integer
            k: MMR TNO integer
            longitude_neptune: heliocentric J2000 longitude of neptune at Epoch

        """
        # initialize the internal variables so they are empty.
        self._a = self._e = self._inc = self._node = self._peri = self._M = None
        self._H = self._lc_gb = self._lc_phase = self._lc_period = self._lc_amplitude = None
        self._phi = self._resamp = self._colors = self._comp = None

        if seed is None:
            seed = numpy.random.randint(1, 999999999)
        self.seed = seed
        self.size = size
        self.distributions = distributions.Distributions(self.seed, self.size)
        self.h_distribution = HDistribution(self.distributions.power_knee_divot,
                                            **dict([('alpha_bright', 1.1),
                                                    ('alpha_faint', 0.4),
                                                    ('h_break', 7.5),
                                                    ('h_max', 11),
                                                    ('h_min', 1)]))
        self.j = j
        self.k = k
        self.longitude_neptune = longitude_neptune  # Neptune's mean longitude on 1 Jan 2013
        self.a_neptune = 30.07 * units.au
        self.comp = "Cls"
        self.epoch = epoch

        # these are protected variables
        self.iter = None  # holds an iterator based on the iterable QTable that holds element values
        self.targets = None  # holds the value of the property target

    @property
    @abstractmethod
    def a(self):
        """
        J2000 Heliocentric semi-major axis.
        """
        pass

    @property
    @abstractmethod
    def e(self):
        """
        J2000 Heliocentric eccentricity
        """
        pass

    @property
    @abstractmethod
    def inc(self):
        """
        J2000 heliocentric inclination of orbit
        """
        pass

    @property
    def node(self):
        """
        Return uniformly distributed nodes.
        """
        if self._node is None:
            self._node = self.distributions.uniform(0, 2 * numpy.pi) * units.rad
        return self._node

    @property
    def peri(self):
        """
        Distribute peri centre to obey the phi/M/_longitude_neptune constraints.
        """
        if self._peri is None:
            self._peri = self.distributions.uniform(0, 2*numpy.pi) * units.rad
        return self._peri

    @property
    def M(self):
        """
        Return uniformly distributed mean anomalies
        """
        if self._M is None:
            self._M = self.distributions.uniform(0, 2*numpy.pi) * units.rad
        return self._M

    @property
    def H(self):
        """A distribution of H values"""
        # define the default size distribution parameters.
        if self._H is None:
            self._H = self.h_distribution()
        return self._H

    @property
    def epoch(self):
        """
        Epoch, in Julian Days, of the orbital elements.
        """
        return self._epoch

    @epoch.setter
    def epoch(self, value):
        if isinstance(value, Quantity):
            value = value.to('day').value
        if isinstance(value, Time):
            value = value.jd
        self._epoch = value * units.day
        # self._epoch = self.distributions.constant(value.to('day').value) * units.day

    @property
    def phi(self):
        """
        Returns the centre of libration for resonances.. 0 otherwise.
        """
        if self._phi is None:
            self._phi = self.distributions.constant(0.0) * units.deg
        return self._phi

    @property
    def resamp(self):
        """
        Returns resamp a 0.0 as this is not a resonant orbit.
        """
        if self._resamp is None:
            self._resamp = self.distributions.constant(0.0) * units.deg
        return self._resamp

    @property
    def j(self):
        """
        Return 0 for non-resonant orbit.
        """
        if self._j is None:
            self._j = self.distributions.constant(0)
        return self._j

    @property
    def k(self):
        """
        Return 0 for non-resonant orbit
        """
        if self._k is None:
            self._k = self.distributions.constant(0)
        return self._k

    @j.setter
    def j(self, value):
        """
        Set the j of j/k resonance description.
        """
        self._j = self.distributions.constant(value)

    @k.setter
    def k(self, value):
        """
        Set the k of j/k resonance description.
        """
        self._k = self.distributions.constant(value)

    @property
    def comp(self):
        """
        Label for the component being generated.
        """
        if self._comp is None:
            self._comp = ['classical', ] * self.size
        return self._comp

    @comp.setter
    def comp(self, value):
        self._comp = [value, ] * self.size

    @property
    def lc_gb(self):
        """
        Opposition surge effect as define in Bowell
        """
        if self._lc_gb is None:
            self._lc_gb = self.distributions.constant(definitions.LIGHT_CURVE_PARAMS['gb']) * units.mag
        return self._lc_gb

    @property
    def lc_phase(self):
        """
        Phase of lightcurve at self.epoch
        """
        if self._lc_phase is None:
            self._lc_phase = self.distributions.constant(definitions.LIGHT_CURVE_PARAMS['phase']) * units.rad
        return self._lc_phase

    @property
    def lc_period(self):
        """
        period of lightcurve
        """
        if self._lc_period is None:
            self._lc_period = self.distributions.constant(definitions.LIGHT_CURVE_PARAMS['period']) * units.day
        return self._lc_period

    @property
    def lc_amplitude(self):
        """
        peak-to-peak amplitude of lightcurve
        """
        if self._lc_amplitude is None:
            self._lc_amplitude = self.distributions.constant(definitions.LIGHT_CURVE_PARAMS['amplitude']) * units.mag
        return self._lc_amplitude

    @property
    def colors(self):
        """
        colors of objects expressed as a list of list:

        [ [ (g-x), (g-r), (i-x), (z-x), (u-x), (V-x), (B-x), (R-x), (I-x), (w-x) ], ... ]
        """
        if self._colors is None:
            g_x = self.distributions.normal(0.7, 0.2) * units.mag  # mu 0.7 and std 0.2 (from CFEPS data)
            r_x = self.distributions.constant(0) * units.mag
            i_x = self.distributions.constant(0) * units.mag
            z_x = self.distributions.constant(0) * units.mag
            u_x = self.distributions.constant(0) * units.mag
            V_x = self.distributions.constant(0) * units.mag
            B_x = self.distributions.constant(0) * units.mag
            R_x = self.distributions.constant(-0.1) * units.mag
            I_x = self.distributions.constant(0) * units.mag
            w_x = self.distributions.constant(0) * units.mag

            self._colors = numpy.column_stack((g_x, r_x, i_x, z_x, u_x, V_x, B_x, R_x, I_x, w_x))
        return self._colors

    def _generate_targets(self):
        """
        Generate the orbit elements and properties of a list of targets to be passed to the simulator.

        Expected to return a QTable or dictionary.  If a QTable then len of table should be self.size.
        If dictionary then each dictionary key should point to a list of length self.size.

        All values stored as Quantity objects to allow conversion to desired units before passing to the SurveySubsF95.detos1

        Must define at least {'a': [], 'e': [], 'inc': [], 'node': [], 'peri': [], 'M': [], 'H': []}  see OSSSSim.simulate for full
        list of keys that can be returned.

        Returns:
            (QTable or dict): set of Quantity objects describing targets.
        """
        self._a = self._e = self._inc = self._node = self._peri = self._M = None
        self._H = self._lc_gb = self._lc_phase = self._lc_period = self._lc_amplitude = None
        self._resamp = self._colors = None

        return QTable([self.a, self.e, self.inc, self.node, self.peri, self.M, [self.epoch,]*len(self.M),
                       self.H, self.lc_gb, self.lc_phase, self.lc_period, self.lc_amplitude,
                       self.resamp, self.colors, self.comp, self.j, self.k, self.phi],
                      names=['a', 'e', 'inc', 'node', 'peri', 'M', 'epoch',
                             'H', 'lc_gb', 'lc_phase', 'lc_period', 'lc_amplitude',
                             'resamp', 'colors', 'comp', 'j', 'k', 'phi'])

    @property
    def targets(self):
        """
        A table of length self.size holding the generated distribution of targets.

        """
        if self._targets is None:
            self._targets = self._generate_targets()
        return self._targets

    @targets.setter
    def targets(self, value):
        """
        set targets to Table stored in value.

        Args:
            value (Table or None): value to set targets to.
        """
        if not isinstance(value, QTable) and value is not None:
            raise ValueError(f"Attempted to set targets table to something that isn't a table: {value}")
        self._targets = value

    @property
    def iter(self):
        """
        An iterator on self.targets table
        """
        if self._iter is None:
            self._iter = iter(self.targets)
        return self._iter

    @iter.setter
    def iter(self, value):
        self._iter = value

    def __iter__(self):
        return self

    def __next__(self):
        """
        return the next row from the orbits table via the `iterator` on the table.

        If we hit the end of orbits table then call draw_distribution to refresh the table.
        """
        try:
            row = next(self.iter)
        except StopIteration:
            # Clear the targets table so a new distribution will be generated.
            self._iter = None
            self.targets = None
            row = next(self.iter)

        return row


class Resonant(Parametric):
    """
    This class defines methods needed to build a parametric Outer Solar System model for use as a model input for OSSSSim
    """

    def __init__(self,
                 size=10**6,
                 seed=123456789,
                 component='Res',
                 longitude_neptune=5.876 * units.rad,
                 epoch=2456839.5 * units.day,
                 j=0,
                 k=0,
                 res_amp_low=0*units.deg,
                 res_amp_mid=5*units.deg,
                 res_amp_high=10*units.deg,
                 res_centre=0*units.deg):
        """
        Setup the boundaries of the simulation.  size and seed are used to initialize a dist_utils.Distribution class.

        longitude_neptune and epoch are stored in self,longitude_neptune and self.epoch for use in
        self._generate_targets.

        See examples/models.py for an example implementation.

        Args:
            size (int): Determines size of the arrays to be generated (default=10^6).
            seed (int): Initialize distributions with this seed to enable reproducible models.
            component (str): Name of this component, stored in output files.
            longitude_neptune (Quantity): The ecliptic longitude of Neptune at epoch_neptune of elements.
            epoch (Quantity): The epoch_neptune of the given longitude of Neptune.
            j (int): the j/k MMR with Neptune, set to 0 if not a resonant orbit.
            k (int): the j/k MMR with Neptune, set to 0 if not a resonant orbit.
            res_centre (Quantity): centre of the resonance phi libration.
            res_amp_low (Quantity): low end of amplitude of phi oscillation
            res_amp_mid (Quantity): middle of the amplitude of phi oscillation
            res_amp_high (Quantity): top end of phi oscillation

            For resonant amplitude is drawn from a distribution that starts at res_amp_low,
            peaks at res_amp_mid and back to 0 at res_amp_high.  Generally uses distributions.Distribution.triangle

        """
        super().__init__(size=size, seed=seed, epoch=epoch, longitude_neptune=longitude_neptune,
                         j=j, k=k)

        if j is None or k is None:
            ValueError(f"Resonance j/k are not None for Resonant Model objects")
        self._res_amp_low = self._res_amp_high = self._res_amp_mid = self._res_centre = None
        self.res_amp_low = res_amp_low
        self.res_amp_high = res_amp_high
        self.res_amp_mid = res_amp_mid
        self.res_centre = res_centre
        self.comp = component.replace(" ", "_")

    @property
    def res_amp_low(self):
        """Low end of resonance amplitude"""
        return self._res_amp_low

    @res_amp_low.setter
    def res_amp_low(self, value):
        self._res_amp_low = value

    @property
    def res_amp_mid(self):
        """Low end of resonance amplitude"""
        return self._res_amp_mid

    @res_amp_mid.setter
    def res_amp_mid(self, value):
        self._res_amp_mid = value

    @property
    def res_amp_high(self):
        """Low end of resonance amplitude"""
        return self._res_amp_high

    @res_amp_high.setter
    def res_amp_high(self, value):
        self._res_amp_high = value

    @property
    def res_centre(self):
        """Low end of resonance amplitude"""
        return self._res_centre

    @res_centre.setter
    def res_centre(self, value):
        self._res_centre = value

    @property
    def a(self):
        """
        J2000 Heliocentric semi-major axis sampled as +/- 0.5 from the resonance semi-major axis value.
        """
        if self._a is None:
            a0 = (self.a_neptune ** (3 / 2) * self.j / self.k) ** (2 / 3)
            a_min = a0 - 0.5 * units.au
            a_max = a0 + 0.5 * units.au
            self._a = self.distributions.uniform(a_min.to('au').value,
                                                 a_max.to('au').value) * units.au
        return self._a


    @property
    def e(self):
        """
        Set the maximum value of 'e' based on the peri-center location of Neptune, minimum value set to 0.02 then randomly sample this
        range of e.
        """
        if self._e is None:
            self._e = self.distributions.uniform(0.05, 0.25)
        return self._e

    @property
    def inc(self):
        """
        Distribute the inclinations based on Brown 2001 functional form.
        """
        if self._inc is None:
            self._inc = self.distributions.truncated_sin_normal(0,
                                                                numpy.deg2rad(11),
                                                                0,
                                                                numpy.deg2rad(40)) * units.rad
        return self._inc

    @property
    def peri(self):
        """
        Distribute peri centre to obey the phi/M/_longitude_neptune constraints.
        """
        """
        From Volk/Chen code.
        phi52 = (libc + 2*resamp*(ran3(seed) - 0.5))
        peri = phi52 - 2.0d0*m + longitude_neptune  - node
        """
        if self._peri is None:
            self._peri = (self.phi / self.k - self.j * self.M / self.k + self.longitude_neptune - self.node) % (360 * units.deg)
        return self._peri

    @property
    def phi(self):
        """
        Compute the phi, libration centre from the resonance centre and sampling the resonance amplitude via sin() weighting.
        """
        if self._phi is None:
            amplitudes = numpy.sin(self.distributions.uniform(0, 2*numpy.pi))
            self._phi = self.res_centre + amplitudes*self.resamp
        return self._phi

    @property
    def resamp(self):
        """
        amplitude of the distribution of resonance centres around libration centre, used in self.phi
        """
        if self._resamp is None:
            self._resamp = self.distributions.triangle(self.res_amp_low.to('rad').value,
                                                       self.res_amp_mid.to('rad').value,
                                                       self.res_amp_high.to('rad').value) * units.rad
        return self._resamp
