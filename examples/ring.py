"""
Create a ring of objects representing theoretical populations of objects in the distant solar system.

When run as a script, uses CFEPS survey characterization, runs until 28 sources have been detected in this ring and plots the result.
"""
from ossssim.models import Parametric
from ossssim import OSSSSim, DetectFile, ModelFile, definitions, plotter
from astropy import units
from astropy.table import Table
import os

class Ring(Parametric):
    """
    Class used to create and store the objects generated and passed by the GiMeObj module into the main Driver.py
    module that executes the survey simulator code.
    """
    def __init__(self, ring_centre, ring_width, model_band='r', **kwargs):
        """Build a ring distribution of given with at a given distance.  Ring is edge-on and circular.

        Args:
            ring_center (units.Quantity): The location of the ring, given as unit quantity
            ring_width (units.Quantity): Width of the ring.
        """
        super().__init__(model_band=model_band, **kwargs)
        self.ring_center = ring_centre
        self.ring_width = ring_width


    @property
    def a_distribution(self):
        """
        Semi-major axis distribution for a narrow ring
        """
        return self.distributions.normal(self.ring_center.to('au').value,
                                         self.ring_width.to('au').value) * units.au



    @property
    def e_distribution(self):
        """
        Eccentricity axis distribution for a narrow ring
        """
        return self.distributions.constant(0.0)

    @property
    def inc_distribution(self):
        """
        Inclination distribution for a narrow ring
        """
        return self.distributions.constant(0.0) * units.rad


def run(model_filename,
        detect_filename,
        characterization_directory, seed, ntrack):
    """
    Using the ParametricModel defined here run the survey simulator and save detected sources to detect_filename

    Args:
        model_filename (str): Name of file to store the model targets.
        detect_filename (str): Name of file to store detected targets.
        characterization_directory (str): Relative or absolute path to directory on disk where the characterization files are organized
        seed (int): random number seed, specifying allows reproducibility
        ntrack (int): < 0 continue for ntrack iterations;
                      > 0 continue until ntrack tracked detections;
                      = 0 continue until input exhausted
    """

    # the default Resonant class arguments setup for a Plutino model....
    model = Ring(45*units.au, 1*units.au, seed=seed, component='Ring', size=1, H_max=9)
    model_file = DetectFile(model_filename,
                            epoch=model.epoch,
                            longitude_neptune=model.longitude_neptune,
                            seed=model.seed)
    detect_file = DetectFile(detect_filename,
                             epoch=model.epoch,
                             longitude_neptune=model.longitude_neptune,
                             seed=model.seed)

    ssim = OSSSSim(characterization_directory=characterization_directory,
                   seed=model.seed)
    n_iter = n_track = n_hits = 0
    for row in model:
        n_iter += 1
        result = ssim.simulate(row, colors=model.colors, model_band=model.model_band, epoch=model.epoch)
        model_file.write_row(result)
        if result['flag'] > 0:
            n_hits += 1
        if result['flag'] == 4:
            n_track += 1
            detect_file.write_row(result)
        if (0 < ntrack <= n_track) or (0 < -ntrack <= n_iter):
            break

    detect_file.write_footer(n_iter=n_iter, n_hits=n_hits, n_track=n_track)
    model_file.write_footer(n_iter=n_iter,
                            n_hits=n_hits,
                            n_track=n_track)


def face_down_plot(model_file: str, detect_file: str) -> None:
    """_
    Plot the detected objects in a face-down plot
    Args:
        detect_file: name of file with the detected sources
    """
    plot = plotter.RosePlot(definitions.Neptune['Longitude'])
    plot.add_model(ModelFile(model_file), mc='k', ms=0.05, alpha=0.1)
    plot.add_model(ModelFile(detect_file), ms=5, mc='g')
    # plot.show()
    plot.savefig('ring.png')

def delete_file_if_exists(filename):
    if os.access(filename, os.F_OK):
        os.remove(filename)


if __name__ == '__main__':
    from ossssim import Characterizations
    seed = 123456789
    n_track = 10
    for i in range(2):
        file_format='ecsv'
        model_filename = f'RingModel_{i}.{file_format}'
        detect_filename = f'RingDetect_{i}.{file_format}'
        delete_file_if_exists(model_filename)
        delete_file_if_exists(detect_filename)
        run(model_filename, detect_filename,
            Characterizations.surveys['CFEPS'],
            seed,
            n_track)
        print(f"iter  : {i}")
        print(f"niter :{len(Table.read(model_filename, format='ascii.ecsv'))}")
        print(f"ntrack:{len(Table.read(detect_filename, format='ascii.ecsv'))}")
    #face_down_plot('RingModel.dat', 'RingDetect.dat')

