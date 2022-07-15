"""
Generate model KBOs orbits from a distribution function and determine if any are detected by a given survey.
"""
import numpy as np
from astropy import units
from astropy.table import Table, QTable

import osssim
from osssim.parametric import Resonant


def run(detect_filename, characterization_directory, seed, ntrack):
    """
    Using the ParametricModel defined here run the survey simulator and save detected sources to detect_filename

    Args:
        detect_filename (str): Name of file to store detected targets.
        characterization_directory (str): Relative or absolute path to directory on disk where the characterization files are organized
        seed (int): random number seed, specifying allows reproducibility
        ntrack (int): < 0 continue for ntrack iterations;
                      > 0 continue until ntracked tracked detections;
                      = 0 continue until input exhasted
    """
    ssim = osssim.OSSSSim(characterization_directory=characterization_directory)

    # the default Resonant class arguments setup for a Plutino model....
    model = Resonant(seed=seed)
    model.comp = "plut"

    detect_file = osssim.DetectFile(detect_filename)
    detect_file.epoch = model.epoch_neptune
    detect_file.lambda_neptune = model.longitude_neptune
    detect_file.colors = osssim.SSimModelFile.COLORS.values()
    detect_file.write_header(seed)

    n_iter = n_track = n_hits = 0
    for row in model:
        n_iter += 1
        result = ssim.simulate(row, seed=model.seed)
        if result['flag'] > 0:
            n_hits += 1
            detect_file.write_row(result)
        if result['flag'] > 2:
            n_track += 1
        if (0 < ntrack < n_track) or (0 < -ntrack < n_iter):
            break

    detect_file.write_footer(n_iter=n_iter, n_hits=n_hits, n_track=n_track)


if __name__ == '__main__':
    run('SimulDetect.dat', '../../../Surveys/CFEPS', 123456789, 28)
