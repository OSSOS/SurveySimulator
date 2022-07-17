"""
Create a ring of objects representing theoretical populations of objects in the distant solar system.
"""
import osssim
from osssim.parametric import NonResonant
from astropy import units


class Ring(NonResonant):
    """
    Class used to create and store the objects generated and passed by the GiMeObj module into the main Driver.py
    module that executes the survey simulator code.
    """
    def __init__(self, ring_centre, ring_width, **kwargs):
        """Build a ring distribution of given with at a given distance.  Ring is edge-on and circular.

        Args:
            ring_center (units.Quantity): The location of the ring, given as unit quantity
            ring_width (units.Quantity): Width of the ring.
        """
        super(Ring).__init__(**kwargs)
        self.ring_center = ring_centre
        self.ring_width = ring_width

    @property
    def a(self):
        return self.distributions.normal(self.ring_center.to('au'),
                                         self.ring_width.to('au')) * units.au

    @property
    def e(self):
        return self.distributions.constant(0.0)

    @property
    def inc(self):
        return self.distributions.constant(0.0) * units.rad


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
    model = Ring(80*units.au, 1*units.au, seed=seed, component='Ring')

    detect_file = osssim.DetectFile(detect_filename)
    detect_file.epoch = model.epoch
    detect_file.lambda_neptune = None
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
