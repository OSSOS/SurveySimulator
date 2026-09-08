"""
Model classes for the Outer Solar System Survey Simulator.
"""
import copy
import logging
import os
import re
import time
from abc import ABC, abstractmethod
from collections import OrderedDict
from collections.abc import Iterable
import numpy
from astropy import units
from astropy.table import QTable, Table
from astropy.time import Time
from astropy.units import Quantity
from numpy import random
from . import definitions
from . import distributions
from . import color
from .color import PhotSpec
from .core import Cartesian

T_ORB_M_UNITS = definitions.T_ORB_M_UNITS
RE_FLOAT = re.compile('[+-]?\d+\.?\d*[de]?\d*')

INITIAL_TABLE_FORMAT = 'ascii.ecsv'
APPEND_TABLE_FORMAT = 'ascii.no_header'
TABLE_COLUMN_DELIMITER = ','


def get_floats_in_str(line):
    """
    Scan a string (line) and return all float like strings as array of floats.

    Convert floats expressed as fortran doubles like '1.0d0' to '1.0e0' before conversion.
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


class OSSSSimFile(ABC):
    """
    Base class for OSSSSim Parametric model and external model input files
    """
    def __init__(self, filename):
        self.filename = filename
        self._table = None
        self._last_row_written = 0
        self._colors = self._longitude_neptune = None
        self.mask_these_if_not_detected = copy.copy(definitions.observables)

    @property
    @abstractmethod
    def table(self) -> QTable:
        """
        Table of orbital elements and measurements of circumstance of observation (if observed)
        """
        pass

    @property
    @abstractmethod
    def column_names(self) -> list[str]:
        """
        return a list of column names to be written or read from file.
        """
        pass

    @property
    @abstractmethod
    def header(self) -> dict:
        """
        Provide a header dictionary with circumstances of the model.
        Should include epoch, seed, longitude_neptune, colors, and component.
        """
        pass

    @property
    def epoch(self) -> Time:
        """
        Epoch of coordinates of orbit read from model file header.
        """
        if not isinstance(self.header['Epoch'], Time):
            self.header['Epoch'] = Time(self.header['Epoch'], format='jd')
        return self.header['Epoch']

    @epoch.setter
    def epoch(self, value: Time) -> None:
        self.header['Epoch'] = value

    @property
    def seed(self) -> int:
        """
        Seed of the model read from the model file header.
        """
        return int(self.header['Seed'])

    @seed.setter
    def seed(self, value) -> None:
        self.header['Seed'] = value

    @property
    def longitude_neptune(self) -> Quantity:
        """
        Longitude of Neptune at epoch
        """
        value = self.header['Longitude_Neptune']
        if not isinstance(value, Quantity):
            value *= units.rad
        return value

    @property
    def colors(self) -> PhotSpec:
        """Returns color array from file header or default if no color array in header."""
        return PhotSpec(self.header['Colors'])

    @colors.setter
    def colors(self, value):
        self.header['Colors'] = value

    @property
    def model_band(self) -> str:
        return self.header['Model_Band']

    @model_band.setter
    def model_band(self, value):
        self.header['Model_Band'] = value

    def write_row(self, row):
        """
        Append a row to the file.
        """
        self.append(row)
        self.write(self.filename, append=True)

    def append(self, this_row) -> None:
        """
        Given a dictionary of row values append the row to the current table

        :param this_row: Dictionary of values to write to row.
        """
        mask = []
        _table_row = []
        for column_name in self.column_names:
            if column_name not in this_row:
                raise ValueError(f"Column name {column_name} not found in {this_row}.")
            masked = (this_row.get('flag', 0) < 1) & (column_name in self.mask_these_if_not_detected)
            mask.append(masked)
            value = this_row[column_name]
            _table_row.append(value)
        self.table.add_row(_table_row, mask=mask)

    def write(self, filename, append=True,
              column_delimiter=TABLE_COLUMN_DELIMITER,
              table_format=APPEND_TABLE_FORMAT):
        """
        Write the table to a file.
        """
        file_already_exists = os.access(filename, os.F_OK)
        if file_already_exists and not append:
            raise IOError(f"File {filename} already exists and not writing in append mode.")
        if not file_already_exists:
            self._last_row_written = 0
            table_format = INITIAL_TABLE_FORMAT
        column_formats = dict([(column_name,
                                definitions.column_format.get(column_name, definitions.column_format['default']))
                               for column_name in self.column_names])
        with open(self.filename, 'a') as f_obj:
            self.table[self._last_row_written:].write(f_obj,
                                                      delimiter=column_delimiter,
                                                      # meta=self.table.meta,
                                                      format=table_format,
                                                      overwrite=False,
                                                      formats=column_formats)
        self._last_row_written = len(self.table)

    def write_footer(self, n_iter, n_hits, n_track):
        """
        Write a footer with the results of the survey simulation

        This is done as a footer instead of header to allow streaming output.
        """
        with open(self.filename, 'a') as f_detect:
            f_detect.write(f'# Total number of objects:   {n_iter:>11d}\n')
            f_detect.write(f'# Number of detections:      {n_hits:>11d}\n')
            f_detect.write(f'# Number of tracked objects: {n_track:>11d}\n')


class ModelFile(Iterable, ABC):
    """This is an abstract class that selects which type of model file object is needed to read a given input file."""

    def __new__(cls, filename, randomize=False):
        with open(filename, 'r') as f_obj:
            if "%ECSV" in f_obj.readline():
                cls = ModelFileEcsv
            else:
                cls = ModelFileOld
        return super().__new__(cls)

    def close(self):
        self._f_obj.close()


class ModelFileOld(ModelFile):
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
        self._model_band = None
        self._longitude_neptune = None
        self.header_lines = []
        self._f_obj = open(self.filename, 'r')
        self.f_loc = 0
        self._targets = None
        self._f = None

    @property
    def epoch(self) -> Time:
        """
        Epoch of coordinates of orbit read from model file header.
        """
        return Time(float(self.header['JD'][0].replace('d', 'e')), format='jd')

    @epoch.setter
    def epoch(self, value: Time):
        self.header['JD'] = [value.jd]

    @property
    def seed(self):
        return int(self.header.get('Seed', [123456789, ])[0])

    @seed.setter
    def seed(self, value):
        self.header['Seed'] = [value]

    @property
    def longitude_neptune(self):
        """
        Longitude of Neptune at epoch
        """
        if self._longitude_neptune is None:
            default_longitude = definitions.Neptune['Longitude'].to('radian').value
            self._longitude_neptune = float(self.header.get('lambdaN', [default_longitude,])[0]
                                            .replace('d', 'e')) * units.radian
        return self._longitude_neptune

    @staticmethod
    def get_key_of_smallest_abs_value_in_dict(d: dict):
        return

    @property
    def colors(self):
        """Returns color array from file header or default if no color array in header."""
        if self._colors is not None:
            return self._colors
        # pickup the default colors to get the correct order of band differences.
        colors_list = numpy.array([get_floats_in_str(color_str)[0]
                                   for color_str in self.header.get('Colors', '')]) * units.mag
        self._colors = color.PhotSpec.from_old_style_list(colors_list)
        return self._colors

    @property
    def model_band(self) -> str:
        if self._model_band is None:
            # set the model band pass to the minimum value in the default color dictionary
            dictionary_of_band_pass_ratio_values = self.colors.colors['default']
            colors = numpy.array([x.to('mag').value for x in dictionary_of_band_pass_ratio_values.values()])
            self._model_band = list(dictionary_of_band_pass_ratio_values)[numpy.argmin(numpy.fabs(colors))].split('-')[1]
        return self._model_band

    @property
    def model_band_pass(self) -> str:
        return self.model_band

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
                        value = line.split('=')[1].strip().split()
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
            if definitions.column_unit.get(colname, None) is not None:
                value = value * definitions.column_unit[colname]
            row[colname] = value
        return row
    @property
    def table(self):
        return self.targets

    @property
    def targets(self):
        """
        targets set by looping over the entire file and returning a 'QTable'.
        This can be used when you want access to the
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


class ModelFileEcsv(ModelFile, OSSSSimFile):
    """
    A class to drive the SSim using a standard model file.

    ModelFile opens file and reads the header for the epoch, seed, longitude_neptune and colors
    and then loops over or randomly offsets into the file to read model objects.
    """

    def __init__(self, filename, randomize=False):
        super().__init__(filename)
        self.filename = filename
        self._table = None
        if randomize:
            DeprecationWarning("Randomize is no longer supported.")
        self.randomize = False
        self._header = None
        self._header_parsed = False
        self._colnames = None
        self._colors = None
        self._epoch = None
        self._longitude_neptune = None
        self.header_lines = []
        self._f_obj = open(self.filename, 'r')
        self.f_loc = 0
        self._targets = None
        self._f = None
        self._last_row_written = 0

    @classmethod
    def read(cls, filename) -> 'ModelFileEcsv':
        model_file = cls(filename)
        model_file._load_table_data()
        return model_file

    def _load_table_data(self) -> None:
        self._table = QTable.read(self.filename, format=INITIAL_TABLE_FORMAT)
        self._last_row_written = len(self._table)

    def __iter__(self):
        return iter(self.table)

    def __len__(self):
        return len(self.table)

    @property
    def header(self) -> dict:
        return self.table.meta

    @property
    def column_names(self) -> list[str]:
        return self.table.colnames

    @property
    def table(self) -> QTable:
        if self._table is None:
            self._load_table_data()
        return self._table

    @property
    def targets(self):
        """
        targets set by looping over the entire file and returning a 'QTable'.
        This can be used when you want access to the
        entire table of data rather than just reading one-line at a time.
        """
        return self.table


class ResultsFile(OSSSSimFile):
    """
    ModelFile structure for output file from Simulator detections.
    """

    @property
    def column_names(self) -> list[str]:
        return self.table.colnames

    @property
    def header(self) -> dict:
        return self.table.meta

    colnames = ['a', 'e', 'inc', 'node', 'peri', 'M', 'H',
                'q', 'r', 'm_rand', 'H_rand', 'x', 'y', 'z',
                'band', 'comp', 'j', 'k']

    def __init__(self, filename: str, seed: int = None, epoch: Time = definitions.Neptune['Epoch'],
                 longitude_neptune: Quantity = definitions.Neptune['Longitude'], colors: PhotSpec = PhotSpec(),
                 model_band: str = definitions.DEFAULT_MODEL_BAND, randomize=False):
        super().__init__(filename)
        self.meta = dict([('Seed', seed),
                          ('Epoch', epoch),
                          ('Longitude_Neptune', longitude_neptune),
                          ('Colors', colors.colors),
                          ('Creation_time', time.strftime("%Y-%m-%dT%H:%M:%S")),
                          ('Model_Band', model_band)])
        if randomize:
            DeprecationWarning("Randomize is no longer supported.")
        self.filename = filename
        self._table = None

    @property
    def table(self) -> QTable:
        if self._table is not None:
            return self._table
        dtypes = []
        column_units = []
        description = {}
        for column_name in self.colnames:
            dtypes.append(definitions.column_dtype.get(column_name,
                                                       definitions.column_dtype['default']))
            column_units.append(definitions.column_unit.get(column_name, units.dimensionless_unscaled))
            description[column_name] = definitions.column_description.get(column_name, None)
        self._table = QTable(names=self.colnames,
                             masked=True,
                             units=column_units,
                             dtype=dtypes,
                             meta=self.meta,
                             descriptions=description)
        return self._table


class ModelOutputFile(ResultsFile):
    """
    Output format used to store the input model, used when model is parametric,
    and you want to keep a record of input for diagnostics
    """

    colnames = ['a', 'e', 'inc', 'node', 'peri', 'M', 'H', 'q',
                'comp', 'j', 'k', 'x', 'y', 'z']

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        for element in ['x', 'y', 'z']:
            self.mask_these_if_not_detected.remove(element)


class DetectFile(ResultsFile):
    """
    Detected object output file structure.
    """
    colnames = ['a', 'e', 'q', 'inc', 'j', 'k', 'node', 'peri', 'M', 'H',  'band', 'color', 'comp',
                'flag', 'Survey', 'eff',
                'm_int', 'm_rand', 'h_rand', 'Mt',
                'RA', 'DEC', 'r', 'delta', 'x', 'y', 'z']


class FakeFile(ResultsFile):
    """
    List of positions of artificial objects to add to the file.
    """
    colnames = ['a', 'e', 'inc', 'node', 'peri', 'M', 'H', 'mag', 'dra', 'ddec', 'RA', 'DEC']


class TrackFile(ResultsFile):
    """
    Tracked object output file structure.
    """
    colnames = ['a', 'e', 'inc', 'node', 'peri', 'M', 'H', 'q', 'r',
                'Mt', 'm_rand', 'H_rand', 'band', 'Survey', 'comp', 'j', 'k']


class Parametric(OSSSSimFile):
    """
    This abstract class defines methods needed to build a parametric
    Outer Solar System model for use as a model input for OSSSSim
    """

    def __init__(self,
                 size: int = definitions.DEFAULT_SIZE,
                 seed: int = definitions.DEFAULT_SEED,
                 epoch: Time = definitions.Neptune['Epoch'],
                 component: str = 'TNO',
                 longitude_neptune: Quantity = definitions.Neptune['Longitude'],
                 H_min: float = definitions.H_MIN,
                 H_max: float = definitions.H_MAX,
                 model_band: str = definitions.DEFAULT_MODEL_BAND,
                 colors: dict = None,
                 **kwargs) -> None:
        """
        Set up the boundaries of the simulation.  size and seed are used to initialize a dist_utils.Distribution class.

        Args:
            size: Determines size of the arrays to be generated (default=10^6).
            seed: Initialize distributions with this seed to enable reproducible models.
            epoch: JD epoch of elements and Neptune longitude
            j: MMR Neptune integer
            k: MMR TNO integer
            longitude_neptune: heliocentric J2000 longitude of neptune at Epoch

        """
        super().__init__(filename=None)
        # initialize the internal variables so they are empty.
        self._sim = self._longitude_neptune = None
        size = kwargs.get('size', size)
        seed = kwargs.get('seed', seed)
        seed = seed is None and numpy.random.randint(1, 999999999) or seed
        epoch = kwargs.get('epoch', epoch)
        component = kwargs.get('component', component).replace(" ", "_")
        longitude_neptune = kwargs.get('longitude_neptune', longitude_neptune)
        H_min = kwargs.get('H_min', H_min)
        H_max = kwargs.get('H_max', H_max)
        model_band = kwargs.get('model_band', model_band)
        colors = kwargs.get('colors', colors)
        colors = colors is None and PhotSpec() or PhotSpec(colors)
        if 'default' not in colors.colors:
            colors.colors['default'] = PhotSpec.COLORS['default']
        self.meta = dict([('Seed', seed),
                          ('Epoch', epoch),
                          ('Longitude_Neptune', longitude_neptune),
                          ('Colors', colors.colors),
                          ('Creation_time', time.strftime("%Y-%m-%dT%H:%M:%S")),
                          ('Component', component),
                          ('Model_Band', model_band)])
        self.size = size
        self.j_distribution = self.k_distribution = [0,]*self.size
        self.comp = component
        self.H_max = H_max
        self.H_min = H_min
        self.model_band = model_band
        self.orbital_elements = ['a', 'e', 'inc', 'node', 'peri', 'M', 'q', 'H', 'j', 'k', 'phi', 'resamp']
        self._a = self._e = self._inc = self._node = self._peri = self._M = self._q = self._H = self._j = self._k = None
        self.distributions = distributions.Distributions(self.seed, self.size)
        self.a_neptune = definitions.Neptune['SemimajorAxis']
        self.rebound_archive = f"Rebound_Archive.bin"
        self.power_knee_divot_params = dict([('alpha_bright', 1.1),
                                             ('alpha_faint', 0.4),
                                             ('h_break', 7.5),
                                             ('h_max', self.H_max),
                                             ('h_min', self.H_min)])
        self._iter = None
        self._targets = None
        self.cartesian = Cartesian(epoch=self.epoch)
        self._init_elements()

    @property
    def a(self):
        if self._a is None:
            self._a = self.a_distribution
        return self._a

    @property
    def e(self):
        if self._e is None:
            self._e = self.e_distribution
        return self._e

    @property
    def inc(self):
        if self._inc is None:
            self._inc = self.inc_distribution
        return self._inc

    @property
    def q(self):
        if self._q is None:
            self._q = self.q_distribution
        return self._q

    @property
    def M(self):
        if self._M is None:
            self._M = self.M_distribution
        return self._M

    @property
    def peri(self):
        if self._peri is None:
            self._peri = self.peri_distribution
        return self._peri

    @property
    def node(self):
        if self._node is None:
            self._node = self.node_distribution
        return self._node

    @property
    def H(self):
        if self._H is None:
            self._H = self.H_distribution
        return self._H

    @property
    def j(self):
        if self._j is None:
            self._j = self.j_distribution
        return self._j

    @property
    def k(self):
        if self._k is None:
            self._k = self.k_distribution
        return self._k

    @property
    def header(self) -> dict:
        return self.meta

    @property
    def column_names(self) -> list[str]:
        return self.table.colnames

    def _init_elements(self):
        """
        Initialize the distributions to trigger generating a new set of model objects
        """
        for element in self.orbital_elements:
            setattr(self, f"_{element}", None)
        self._iter = self._targets = self._cartesian = self._sim = None

    @abstractmethod
    def a_distribution(self) -> Quantity:
        """
        Semi-major axis of the orbit.
        """
        pass

    @property
    @abstractmethod
    def e_distribution(self) -> Quantity:
        """
        Semi-major axis of the orbit.
        """
        pass

    @property
    @abstractmethod
    def inc_distribution(self) -> numpy.ndarray:
        """
        Semi-major axis of the orbit.
        """
        pass

    @property
    def node_distribution(self) -> Quantity:
        """
        Return uniformly distributed nodes.
        """
        return self.distributions.uniform(0, 2 * numpy.pi) * units.rad

    @property
    def peri_distribution(self) -> Quantity:
        """
        Distribute peri uniformly.
        """
        return self.distributions.uniform(0, 2 * numpy.pi) * units.rad

    @property
    def M_distribution(self) -> Quantity:
        """
        Return uniformly distributed mean anomalies
        """
        return self.distributions.uniform(0, 2 * numpy.pi) * units.rad

    @property
    def H_distribution(self) -> Quantity:
        """A distribution of H values"""
        return self.distributions.power_knee_divot(**self.power_knee_divot_params) * units.mag

    @property
    def q_distribution(self) -> Quantity:
        return self.a*(1-self.e)

    @property
    def comp(self) -> list:
        """
        Label for the component being generated.
        """
        return [self._comp,] * self.size

    @comp.setter
    def comp(self, value: (str, list)):
        self._comp = value

    @property
    def lc_gb(self) -> Quantity:
        """
        Opposition surge effect as define in Bowell
        """
        return self.distributions.constant(definitions.LIGHT_CURVE_PARAMS['gb'].value) * units.mag

    @property
    def lc_phase(self) -> Quantity:
        """
        Phase of lightcurve at self.epoch
        """
        return self.distributions.uniform(0, 2*numpy.pi) * units.rad

    @property
    def lc_period(self) -> Quantity:
        """
        period of lightcurve
        """
        return self.distributions.uniform(0, definitions.LIGHT_CURVE_PARAMS['period'].value) * units.day

    @property
    def lc_amplitude(self) -> Quantity:
        """
        peak-to-peak amplitude of lightcurve
        """
        return self.distributions.uniform(0, definitions.LIGHT_CURVE_PARAMS['amplitude'].value) * units.mag

    def _generate_targets(self) -> QTable:
        """
        Generate the orbit elements and properties of a list of targets to be passed to the simulator.

        Expected to return a QTable or dictionary.  If a QTable then len of table should be self.size.
        If dictionary then each dictionary key should point to a list of length self.size.

        All values stored as Quantity objects to allow conversion to desired units
        before passing to the ossssimlib.surveysub.detos1

        Must define at least {'a': [], 'e': [], 'inc': [], 'node': [], 'peri': [], 'M': [], 'H': []}
        see sim for full list of keys that can be returned.

        Returns:
            (QTable or dict): set of Quantity objects describing targets.
        """
        rows = {}
        for element in ['a', 'e', 'inc', 'node', 'peri', 'M', 'q', 'H', 'j', 'k', 'comp',
                        'lc_gb', 'lc_phase', 'lc_period', 'lc_amplitude']:
            rows[element] = getattr(self, element)

        pos_zeros = numpy.zeros(self.size)*units.au
        vel_zeros = numpy.zeros(self.size)*units.au/units.yr
        carts = ['x', 'y', 'z']
        for cart in carts:
            rows[cart] = pos_zeros
        carts = ['vx', 'vy', 'vz']
        for cart in carts:
            rows[cart] = vel_zeros
        table = QTable(self.cartesian(rows=rows))
        return table

    @property
    def table(self):
        return self.targets

    @property
    def targets(self):
        """
        A table of length self.size holding the generated distribution of targets.

        """
        if self._targets is None:
            self._targets = self._generate_targets()
        return self._targets

    @property
    def iter(self):
        """
        An iterator on targets table
        """
        if self._iter is None:
            self._iter = iter(self.targets)
        return self._iter

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
            self._init_elements()
            row = next(self.iter)
        return dict(row)


class Implanted(Parametric, ABC):
    """"
    Objects that were implanted into the Kuiper belt region and appear to share a SFD.
    """
    def __init__(self, sigma_i=12, **kwargs):
        super().__init__(**kwargs)
        self.sigma_i = sigma_i

    @property
    def h_distribution(self):
        return distributions.HDistribution(self.distributions.implanted_sfd,
                                           **dict([('h_max', self.H_max),
                                                   ('h_min', self.H_min)]))

    @property
    def inc_distribution(self):
        """
        Distribute the inclinations based on Brown 2001 functional form.
        """
        return self.distributions.truncated_sin_normal(0,
                                                       numpy.deg2rad(self.sigma_i),
                                                       0,
                                                       numpy.deg2rad(45)) * units.rad


class Resonant(Implanted):
    """
    This class defines methods needed to build a parametric
    Outer Solar System model for use as a model input for OSSSSim
    """

    def __init__(self,
                 size=10 ** 5,
                 seed=123456789,
                 comp='Res',
                 j=None,
                 k=None,
                 res_amp_low=0 * units.deg,
                 res_amp_mid=5 * units.deg,
                 res_amp_high=10 * units.deg,
                 res_centre=180 * units.deg, **kwargs):
        """
        Set up the boundaries of the simulation.  size and seed are used to initialize a dist_utils.Distribution class.

        longitude_neptune and epoch are stored in self,longitude_neptune and self.epoch for use in
        self._generate_targets.

        See examples/models.py for an example implementation.


                res_amp_low=20*units.degree,
                res_amp_mid=95*units.degree,
                res_amp_high=130*units.degree,

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
        super().__init__(size=size,
                         seed=seed,
                         comp=comp, **kwargs)
        if j is None or k is None:
            ValueError(f"Resonance j/k cannot be None for Resonant Model objects")
        self.j_distribution = self.distributions.constant(j)
        self.k_distribution = self.distributions.constant(k)
        self._phi = self._resamp = None
        self.orbital_elements.extend(['resamp', 'phi', 'j', 'k'])
        self.res_amp_low = res_amp_low
        self.res_amp_high = res_amp_high
        self.res_amp_mid = res_amp_mid
        self.res_centre = res_centre

    @property
    def res_amp_low(self) -> Quantity:
        """Low end of resonance amplitude"""
        return self._res_amp_low

    @res_amp_low.setter
    def res_amp_low(self, value: Quantity):
        self._res_amp_low = value

    @property
    def res_amp_mid(self) -> Quantity:
        """Low end of resonance amplitude"""
        return self._res_amp_mid

    @res_amp_mid.setter
    def res_amp_mid(self, value: Quantity):
        self._res_amp_mid = value

    @property
    def res_amp_high(self) -> Quantity:
        """Low end of resonance amplitude"""
        return self._res_amp_high

    @res_amp_high.setter
    def res_amp_high(self, value: Quantity):
        self._res_amp_high = value

    @property
    def res_centre(self) -> Quantity:
        """Centre of resonance, aka phi0"""
        return self._phi0

    @res_centre.setter
    def res_centre(self, value: Quantity):
        self._phi0 = value

    @property
    def a_distribution(self) -> Quantity:
        """
        J2000 Heliocentric semi-major axis sampled as +/- 0.5 from the resonance semi-major axis value.
        """
        a0 = (self.a_neptune ** (3 / 2) * self.j[0] / self.k[0]) ** (2 / 3)
        a_min = a0 - 0.5 * units.au
        a_max = a0 + 0.5 * units.au
        return self.distributions.uniform(a_min.to('au').value,
                                          a_max.to('au').value) * units.au

    @property
    def e_distribution(self) -> numpy.array:
        """
        Set the maximum value of 'e' based on the peri-center location of Neptune,
        minimum value set to 0.02 then randomly sample this
        range of e.
        """
        return self.distributions.uniform(0.19, 0.2)

    @property
    def resamp(self):
        """
        amplitude of the distribution of resonance centres around libration centre, used in self.phi
        """
        if self._resamp is None:
            self._resamp = self.resamp_distribution
        return self._resamp

    @property
    def phi(self):
        """
        Libration angle of the resonance.
        """
        if self._phi is None:
            self._phi = self.phi_distribution
        return self._phi

    @property
    def M_distribution(self) -> Quantity:
        """
        Return uniformly distributed mean anomalies
        """
        return self.distributions.uniform(0, self.k[1]*2*numpy.pi) * units.rad

    @property
    def peri_distribution(self) -> Quantity:
        """
        Distribute peri centre to obey the phi/M/_longitude_neptune constraints.

        See Volk et al. 2016 for info on computing peri given a choice phi.
        """
        # below is different algebra to get the same result
        # self._peri = (self.phi / self.k - self.j * self.M / self.k + self.longitude_neptune - self.node)
        # % (360 * units.deg)
        # self._peri = (self.phi - p*self.M + q*self.longitude_neptune - q*self.node)/q
        # return (self.phi - self.j * self.M + self.k * (self.longitude_neptune - self.node)) / self.k
        return ((self.phi - self.j * self.M)/self.k + self.longitude_neptune - self.node) % (360 * units.deg)

    @property
    def phi_distribution(self) -> Quantity:
        """
        Compute the phi, libration centre from the resonance centre and sampling the
        resonance amplitude via sin() weighting.
        """
        amplitudes = numpy.sin(self.distributions.uniform(0, 2 * numpy.pi))
        return self.res_centre + amplitudes * self.resamp

    @property
    def resamp_distribution(self) -> Quantity:
        """
        amplitude of the distribution of resonance centres around libration centre, used in self.phi
        """
        return self.distributions.triangle(self.res_amp_low.to('rad').value,
                                           self.res_amp_mid.to('rad').value,
                                           self.res_amp_high.to('rad').value) * units.rad
