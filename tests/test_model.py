"""
Create a ring of objects representing theoretical populations of objects in the distant solar system.

When run as a script, uses CFEPS survey characterization, runs until 28 sources have been detected in this ring and plots the result.
"""
import ossssim
from ossssim.models import Parametric, Resonant
from ossssim import OSSSSim, DetectFile, PhotSpec, Characterizations
from astropy import units
import unittest
from tempfile import NamedTemporaryFile
import pathlib

script_directory = pathlib.Path(__file__).parent.resolve()


class Plutino(Resonant):
    def __init__(self, **kwargs):
        super().__init__(j=3, k=2, comp="plutino", **kwargs)


class Ring(Parametric):
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
        super().__init__(**kwargs)
        self.ring_center = ring_centre
        self.ring_width = ring_width

    @property
    def a_distribution(self):
        """
        Semi-major axis distribution for a narrow ring
        """
        return self.distributions.constant(40.0) * units.au

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


class ConfirmResonantTest(unittest.TestCase):

    def setUp(self):
        self.characterization_directory = Characterizations.surveys['CFEPS']
        print(self.characterization_directory)
        self.model = Plutino(size=10000)
        self.sim = OSSSSim(characterization_directory=self.characterization_directory, seed=123456789)

    def test_semimajor_axis_limits(self):
        self.assertLess(self.model.a.max(), 39.95*units.au)
        self.assertGreater(self.model.a.min(), 38.90*units.au)


class CreateModelFileTest(unittest.TestCase):
    def setUp(self):
        self.characterization_directory = Characterizations.surveys['CFEPS']
        self.ssim = OSSSSim(characterization_directory=self.characterization_directory,  seed=123456789)
        self.model = Ring(45*units.au, 1*units.au, component='Ring', size=1, H_max=9)
        self.model_filename = NamedTemporaryFile().name

    def test_run(self):
        model_file = DetectFile(self.model_filename,
                                longitude_neptune=self.model.longitude_neptune,
                                epoch=self.model.epoch,
                                colors=PhotSpec())
        for row in self.model:
            self.assertAlmostEqual(row['a'].value, 40.0)
            self.assertAlmostEqual(row['e'], 0.0)
            self.assertAlmostEqual(row['inc'].value, 0.0)
            result = self.ssim.simulate(row,
                                        colors=model_file.colors,
                                        model_band=model_file.model_band,
                                        seed=self.model.seed,
                                        epoch=self.model.epoch)
            model_file.write_row(result)
            if result['flag'] > 0:
                break

        model_file.write_footer(n_iter=1,
                                n_hits=0,
                                n_track=0)


if __name__ == '__main__':
    unittest.main()
