"""
Read in an OSSOS SurveySimulator formated file into an astropy Table.
"""
import logging
import os
import re
import time
from collections import OrderedDict, Iterable

from astropy import units as u
from astropy.time import Time
from astropy.units import Quantity
from numpy.random import random

RE_FLOAT = re.compile('\d+\.?\d*[de]?\d*')

colunits = {'a': u.au,
            'e': None,
            'inc': u.degree,
            'q': u.au,
            'r': u.au,
            'm': u.degree,
            'node': u.degree,
            'peri': u.degree,
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


class SSimResults:
    """
    File structure for output file from Simulator detections.
    """
    COLUMN_WIDTH = 10
    COLUMN_MAP = {'i': 'inc'}

    formats = {'flag': f'{COLUMN_WIDTH}d',
               'Survey': f'{COLUMN_WIDTH}s',
               'Comments': f'{COLUMN_WIDTH}s',
               'comp': f'{COLUMN_WIDTH}s',
               'RA': f'{COLUMN_WIDTH}.5f',
               'DEC': f'{COLUMN_WIDTH}.4f',
               'j': f'{COLUMN_WIDTH}d',
               'k': f'{COLUMN_WIDTH}d',
               'default': f'{COLUMN_WIDTH}.4f'}

    colnames = ['a', 'e', 'inc', 'node', 'peri', 'm', 'h',
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
        self._lambda_neptune = None
        self._seed = None
        self.header_lines = []
        self.floc = None

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
        if not isinstance(value, float):
            raise ValueError(f"setting epoch to non flaot value, expected JD as float")
        self._epoch = Time(value)

    @property
    def lambda_neptune(self):
        """
        Longitude of Neptune at epoch
        """
        if self._lambda_neptune is None:
            raise ValueError(f"lambda_neptune has not been set.")
        return self._lambda_neptune

    @lambda_neptune.setter
    def lambda_neptune(self, value):
        if not isinstance(value, Quantity):
            raise ValueError(f"lambda_neptune must be set to a Quantity with units")
        self._lambda_neptune = value

    @property
    def colors(self):
        """Returns color array from file header or default if no color array in header."""
        if self._colors is None:
            raise ValueError(f"colors list has not been initialized.")
        return self._colors

    @colors.setter
    def colors(self, values):
        if not isinstance(values, list):
            raise ValueError(f"Attempt to set color list with non-list object.")
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
                col_format = self.formats.get(colname, self.formats['default'])
                if isinstance(this_row[colname], Quantity):
                    value = this_row[colname].to(colunits[colname]).value
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
            f_detect.write(f"# Longitude of Neptune: lambdaN = {self.lambda_neptune}\n")
            color_str = [f"{c:5.2f} " for c in self.colors]
            f_detect.write(f"Colors = {color_str}\n")
            f_detect.write(f"#\n")
            dati = time.strftime("%Y-%m-%dT%H:%M:%S.000  %z")
            f_detect.write(f"# Creation_time: {dati:30s}\n")
            f_detect.write('# flag: >0: detected; >2: characterized; 0 mod(2): tracked\n')
            f_detect.write('# Survey: name of the observing block where model object was detected\n')
            f_detect.write("# ")
            for colname in self.colnames:
                f_detect.write("{colname:>{width}s} ".format(colname=colname,
                                                             width=self.COLUMN_WIDTH))
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


class DetectFile(SSimResults):
    """
    Detected object output file structure.
    """
    colnames = ['a', 'e', 'inc', 'node', 'peri', 'm', 'h', 'q', 'r', 'm', 'm_rand', 'h_rand', 'color', 'flag',
                'delta', 'm_int', 'eff', 'RA', 'DEC', 'Survey', 'comp', 'j', 'k']


class TrackFile(SSimResults):
    """
    Tracked object output file structure.
    """
    colnames = ['a', 'e', 'inc', 'node', 'peri', 'm', 'h', 'q', 'r',
                'm_rand', 'H_rand', 'color', 'Survey', 'comp', 'j', 'k']


class SSimModelFile(Iterable):
    """
    A class to drive the SSim using a standard model file.

    SSimModelFile opens file and reads the header for the epoch, seed, lamda_neptune and colors
    and then loops over or randomly offsets into the file to read model objects.
    """

    # Default color values, appropriate for mean values of TNOs
    COLORS = OrderedDict(
        (('g-x', +0.0),
         ('r-x', -0.7),
         ('i-x', -1.2),
         ('z-x', -1.7),
         ('u-x', +0.8),
         ('V-x', +0.5),
         ('B-x', +0.1),
         ('R-x', -0.8),
         ('I-x', -1.2),
         ('x-x', +0.0)))

    COLUMN_MAP = {'i': 'inc'}

    def __init__(self, filename, randomize=False):
        self.filename = filename
        self.randomize = randomize
        self._header = None
        self._header_parsed = False
        self._colnames = None
        self._colors = None
        self._epoch = None
        self._lambda_neptune = None
        self._seed = None
        self.header_lines = []
        self.floc = None
        self._fobj = open(self.filename, 'r')

    @property
    def epoch(self):
        """
        Epoch of coordinates of orbit read from model file header.
        """
        if self._epoch is None:
            self._epoch = Time(float(self.header['JD'].replace('d', 'e')),
                               format='jd')
        return self._epoch

    @property
    def lambda_neptune(self):
        """
        Longitude of Neptune at epoch
        """
        if self._lambda_neptune is None:
            self._lambda_neptune = float(self.header['lambdaN'].replace('d', 'e')) * u.radian
        return self._lambda_neptune

    @property
    def colors(self):
        """Returns color array from file header or default if no color array in header."""
        if self._colors is None:
            # pickup the
            self._colors = list(self.COLORS.values())
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
                colname = self.COLUMN_MAP.get(colname, colname)
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
        with open(self.filename, 'r') as fobj:
            if self.floc is not None:
                fobj.seek(self.floc)
            self.floc = fobj.tell()
            for line in fobj.readlines():
                # line = line.decode('utf-8')
                if line.startswith('#'):
                    logging.debug(f"Parsing Comment: {line}")
                    self.header_lines.append(line[1:])
                    if line.strip() == "#":
                        continue
                    if '=' in line:
                        keyword = line.split('=')[0].split()[-1]
                        value = line.split('=')[1].strip()
                        self._header[keyword] = value
                    previous_line = line
                    self.floc = fobj.tell()
                else:
                    if previous_line is not None:
                        # expect that the last header line is the column name header.
                        self._header['colnames'] = previous_line[1:]
                    break
            self._header_parsed = True
        if self._header is None:
            raise IOError(f"Failed to parse keywods from header of {self.filename}")
        self._header['_end_of_header_offset'] = self.floc
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
                self._fobj.seek(random.randint(self.header['_end_of_header_offset'],
                                               os.stat(self.filename).st_size))
                try:
                    # read to the end of this line.
                    self._fobj.readline()
                    while True:
                        line = self._fobj.readline()
                        if len(line) == 0:
                            raise EOFError
                        if line[0] != "#":
                            break
                    break
                except EOFError:
                    self._fobj.close()
                    self._fobj = open(self.filename)
                    pass
        else:
            while True:
                line = self._fobj.readline()
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
            if colunits[colname] is not None:
                value = value * colunits[colname]
            row[colname] = value
        return row


def run(filename):
    """
    return the first row of a model file as a test of the SSimModelFile class.
    """
    reader = SSimModelFile(filename)
    for row in reader:
        print(row)
        break


if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG)
    run('../../../Models/L7model-3.0-9.0')
