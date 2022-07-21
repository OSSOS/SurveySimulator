"""
Generate model KBOs orbits from a distribution function and determine if any are detected by a given survey.
"""
from ossssim.models import Resonant, DetectFile
from ossssim import definitions, plotter, ModelFile
from ossssim import OSSSSim


def run(model_filename,
        detect_filename, characterization_directory, seed, n_track):
    """
    Using the ParametricModel defined here run the survey simulator and save detected sources to detect_filename

    Args:
        model_filename (str): Name of file to store model objects.
        detect_filename (str): Name of file to store detected targets.
        characterization_directory (str): Relative or absolute path to directory on disk where the characterization files are organized
        seed (int): random number seed, specifying allows reproducibility
        n_track (int): < 0 continue for ntrack iterations;
                      > 0 continue until n_tracked tracked detections;
                      = 0 continue until input exhausted
    """
    ssim = OSSSSim(characterization_directory=characterization_directory)

    # the default Resonant class arguments setup for a Plutino model....
    model = Resonant(j=3, k=2, seed=seed)
    model.comp = "Plutino"

    model_file = DetectFile(model_filename)
    model_file.epoch = model.epoch
    model_file.lambda_neptune =model.longitude_neptune
    model_file.colors = definitions.COLORS.values()
    model_file.write_header(seed)

    detect_file = DetectFile(detect_filename)
    detect_file.epoch = model.epoch
    detect_file.lambda_neptune = model.longitude_neptune
    detect_file.colors = definitions.COLORS.values()
    detect_file.write_header(seed)

    n_iter = n_tracked = n_hits = 0
    for row in model:
        n_iter += 1
        result = ssim.simulate(row, seed=model.seed)
        model_file.write_row(result)
        if result['flag'] > 0:
            n_hits += 1
            detect_file.write_row(result)
        if result['flag'] > 2:
            n_tracked += 1
        if (0 < n_track < n_tracked) or (0 < -n_track < n_iter):
            break

    detect_file.write_footer(n_iter=n_iter, n_hits=n_hits, n_track=n_tracked)
    model_file.write_footer(n_iter=n_iter, n_hits=n_hits, n_track=n_tracked)

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
    run('PlutinoModel.dat',
        'PlutinoDetect.dat',
        'Surveys/CFEPS', 123456789, 28)

    face_down_plot('PlutinoModel.dat',
                   'PlutinoDetect.dat')
