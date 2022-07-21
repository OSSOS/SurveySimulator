"""
An example of a 2:1 model for use with OSSSSIM.simulate, see also plutino.py for more details
"""
from astropy import units
import numpy
import ossssim
from ossssim.models import Resonant
from ossssim import definitions, plotter, ModelFile


def run(model_filename, detect_filename, characterization_directory, seed, ntrack):
    """
    Using the Resonant defined in ossssim.parametric build model objects and pass them through the survey simulator, save detections

    Args:
        detect_filename (str): Name of file to store detected targets.
        characterization_directory (str): Relative or absolute path to directory on disk where the characterization files are organized
        seed (int): random number seed, specifying allows reproducibility
        ntrack (int): < 0 continue for ntrack iterations;
                      > 0 continue until n_tracked tracked detections;
                      = 0 continue until input exhausted
    """
    ssim = ossssim.OSSSSim(characterization_directory=characterization_directory)

    sym_model = Resonant(j=2, k=1, res_centre=0 * units.deg, component='Symmetric')
    leading_model = Resonant(j=2, k=1, res_centre=80 * units.deg, component='Asymmetric_Leading')
    trailing_model = Resonant(j=2, k=1, res_centre=280 * units.deg, component='Asymmetric_Trailing')

    model_file = ossssim.DetectFile(model_filename)
    model_file.epoch = sym_model.epoch
    model_file.lambda_neptune = sym_model.longitude_neptune
    model_file.colors = definitions.COLORS.values()
    model_file.write_header(seed)

    detect_file = ossssim.DetectFile(detect_filename)
    detect_file.epoch = sym_model.epoch
    detect_file.lambda_neptune = sym_model.longitude_neptune
    detect_file.colors = definitions.COLORS.values()
    detect_file.write_header(seed)

    n_iter = n_track = n_hits = 0
    while True:
        # pick which model to pass, based on symmetric,trailing/leading ratios
        selection = numpy.random.rand()
        if selection < 1/3:
            model = sym_model
        elif selection < 2/3:
            model = leading_model
        else:
            model = trailing_model
        row = next(model)
        n_iter += 1
        result = ssim.simulate(row, seed=model.seed)
        model_file.write_row(result)
        if result['flag'] > 0:
            n_hits += 1
            detect_file.write_row(result)
        if result['flag'] > 2:
            n_track += 1
        if (0 < ntrack < n_track) or (0 < -ntrack < n_iter):
            break

    detect_file.write_footer(n_iter=n_iter, n_hits=n_hits, n_track=n_track)
    model_file.write_footer(n_iter=n_iter, n_hits=n_hits, n_track=n_track)


def face_down_plot(model_file: str, detect_file: str) -> None:
    """_
    Plot the detected objects in a face-down plot
    Args:
        detect_file: name of file with the detected sources
    """
    plot = plotter.FaceDownPlot(definitions.Neptune['Longitude'])
    plot.add_model(ModelFile(model_file), mc='k', ms=0.05, alpha=0.1)
    plot.add_model(ModelFile(detect_file), ms=5, mc='g')
    plot.plot_rings()
    plot.show()


if __name__ == '__main__':
    run('TwotinoModel.dat', 'TwotinoDetect.dat', 'Surveys/CFEPS', 123456789, 28)

    face_down_plot('TwotinoModel.dat',
                   'TwotinoDetect.dat')
