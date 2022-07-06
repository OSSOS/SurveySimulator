"""
Read in an OSSOS SurveySimulator formated file into an astropy Table.
"""
import logging
import re
import time
from astropy import units as u
from astropy.units import Quantity

RE_FLOAT = re.compile('\d+\.?\d*[de]?\d*')


def get_floats_in_str(line):
    """
    Scan a string (line) and return all float like strings as array of floats.

    Converts floats expressed as fortran doubles like '1.0d0' to '1.0e0' before conversion.
    """
    values = RE_FLOAT.findall(line)

    try:
        result = [float(x.replace('d', 'e')) for x in values]
    except ValueError as ex:
        print(line)
        print(values)
        raise ex
    return result

COLUMN_WIDTH = 10

column_units = {
    'a': u.au,
    'e': None,
    'inc': u.degree,
    'q': u.au,
    'node': u.degree,
    'peri': u.degree,
    'm': u.degree,
    'h': u.mag,
    'dist': u.au,
    'comp': None,
    'j': None,
    'k': None
}


class SSimOutFile:
    """
    File structure for output file from Simulator detections.
    """
    formats = {'flag': f'{COLUMN_WIDTH}d',
               'Survey': f'{COLUMN_WIDTH}s',
               'Comments': f'{COLUMN_WIDTH}s',
               'comp': f'{COLUMN_WIDTH}s',
               'RA': f'{COLUMN_WIDTH}.5f',
               'DEC': f'{COLUMN_WIDTH}.4f',
               'j': f'{COLUMN_WIDTH}d',
               'k': f'{COLUMN_WIDTH}d',
               'default': f'{COLUMN_WIDTH}.3f'}
    colnames = ['a', 'e', 'inc', 'q', 'r', 'M', 'node', 'peri', 'm_rand', 'H_rand', 'color', 'comp', 'j', 'k']

    def __init__(self, filename):
        self.filename = filename

    def write_row(self, row):
        """
        Given a dictionary of row values write the row according to the order in colnames

        :param row: Dictionary of values to write to row.
        """
        with open(self.filename, 'a') as f_detect:
            sep = " "
            for colname in self.colnames:
                col_format = self.formats.get(colname, self.formats['default'])
                if isinstance(row[colname], Quantity):
                    value = row[colname].to(DetectFile.colunits[colname]).value
                else:
                    value = row[colname]
                o_str = "{sep}{value:{col_format}}".format(value=value, col_format=col_format, sep=sep)
                f_detect.write(o_str)
                sep = " "
            f_detect.write('\n')

    def write_header(self, seed):
        """
        :param seed: the seed for random generator that resulted in these detections.
        :type seed: int
        """
        with open(self.filename, 'w') as f_detect:
            f_detect.write(f"# Seed: {seed:10d}\n")
            dati = time.strftime("%Y-%m-%dT%H:%M:%S.000  %z")
            f_detect.write(f"# Creation time: {dati:30s}\n")
            f_detect.write('# flag: >0: detected; >2: characterized; 0 mod(2): tracked\n')
            f_detect.write('# Survey: name of the observing block where model object was detected\n')
            sep = "#"
            for colname in self.colnames:
                f_detect.write("{sep}{colname:>{width}s}".format(sep=sep, colname=colname, width=COLUMN_WIDTH))
                sep = " "
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


class DetectFile(SSimOutFile):
    """
    Detected object output file structure.
    """
    colnames = ['a', 'e', 'inc', 'q', 'r', 'm', 'node', 'peri', 'm_rand', 'h_rand', 'color', 'flag',
                     'delta', 'm_int', 'h', 'eff', 'RA', 'DEC', 'Survey', 'comp', 'j', 'k']
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
                'eff': None,
                'RA': u.deg,
                'DEC': u.deg,
                'Survey': None,
                'comp': None,
                'j': None,
                'k': None
                }


class TrackFile(SSimOutFile):
    """
    Tracked object output file structure.
    """
    colnames = ['a', 'e', 'inc', 'q', 'r', 'M', 'node', 'peri',
                'm_rand', 'H_rand', 'color', 'Survey', 'comp', 'j', 'k']




if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG)
    reader = SSimModelReader('../../../Models/L7model-3.0-9.0')
    for row in reader:
        print(row)
        break

