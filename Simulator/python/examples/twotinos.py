"""
An example of a 2:1 model for use with OSSSSIM.simulate, see also toy.py for more details
"""
from astropy import units
import osssim
from osssim.parametric import Resonant
import numpy
from osssim import definitions


def run(detect_filename, characterization_directory, seed, ntrack):
    """
    Using the Resonant defined in osssim.parametric build model objects and pass them through the survey simulator, save detections

    Args:
        detect_filename (str): Name of file to store detected targets.
        characterization_directory (str): Relative or absolute path to directory on disk where the characterization files are organized
        seed (int): random number seed, specifying allows reproducibility
        ntrack (int): < 0 continue for ntrack iterations;
                      > 0 continue until ntracked tracked detections;
                      = 0 continue until input exhasted
    """
    ssim = osssim.OSSSSim(characterization_directory=characterization_directory)

    symm_model = Resonant(j=2, k=1, res_centre=0 * units.deg, component='Symmetric')
    leading_model = Resonant(j=2, k=1, res_centre=80 * units.deg, component='Asym_Leading')
    trailing_model = Resonant(j=2, k=1, res_centre=280 * units.deg, component='Asym_Trailing')

    detect_file = osssim.DetectFile(detect_filename)
    detect_file.epoch = symm_model.epoch_neptune
    detect_file.lambda_neptune = symm_model.longitude_neptune
    detect_file.colors = definitions.COLORS.values()
    detect_file.write_header(seed)

    n_iter = n_track = n_hits = 0
    while True:
        # pick which model to pass, based on symmetric,trailing/leading ratios
        selection = numpy.random.rand()
        if selection < 1/3:
            model = symm_model
        elif selection < 2/3:
            model = leading_model
        else:
            model = trailing_model
        row = next(model)
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
