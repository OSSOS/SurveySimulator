"""
An example of a 2:1 model for use with OSSSSIM.simulate.  This is the OSSOS++ 2:1 model from Chen and Volk.

To work in the python mode here we define two new classes derived from the Resonance.  The Resonance base class handles computing
peri from the libration centre and by default will distributed the libration centre using a 'low, mid, high' set of values that define
a triangle distribution (starts at 0 at low, rises to 1 at mid, back to 0 at high) for the amplitude distribution.  That was CFEPS.

 For OSSOS we have a more complex set computations where the libration centre (phi0) and the libration amplitude (resamp) about that
 centre depend on the eccentricity of orbit and type (Symmetric or Asymmetric) of the resonance.
"""
import numpy
from astropy import units
from astropy.coordinates import SkyCoord
from astropy.table import Table
import ossssim
from ossssim.models import Resonant
import os

SURVEYS_DIR = '/arc/projects/classy/Surveys'

class TwoTinoSymmetric(Resonant):
    """
    Build a Symmetric class from the resonant base class, j/k to 2/1.
    """
    def __init__(self, e_min=0.07, e_max=0.30, j=2, k=1, **kwargs):
        """
        Symmetric, libration centers are all 180
        pick libration amplitudes uniformly from 125-165
        pick e uniformly from 0.05 -> 0.3 (Chiang&Jordan02)

        """
        super().__init__(j=j, k=k, **kwargs)
        self.e_min = e_min
        self.e_max = e_max
        self.phi0 = 180.0 * units.deg

    @property
    def e(self):
        """
        Set the maximum value of 'e' based on the peri-center location of Neptune, minimum value set to 0.02 then randomly sample this
        range of e.
        """
        if self._e is None:
            # max is set by q but also limited by users choice of e_max.
            e_max = 1 - ((22 * units.au) / self.a).value
            e_max[e_max > self.e_max] = self.e_max
            self._e = self.e_min + self.distributions.uniform(0, 1) * (e_max - self.e_min)
        return self._e

    @property
    def resamp(self):
        """
        Uniformly distributed across the width
        """
        if self._resamp is None:
            self._resamp = self.distributions.uniform(140., 167.) * units.deg
        return self._resamp

    @property
    def phi(self):
        """
        Resonance centre is the libration centre +/- the resamp

        The default phi is distributed using a sin weighted distribution but in the OSSOS++ model this was switched to uniform for
        the symmetric resonance.
        """
        if self._phi is None:
            self._phi = self.phi0 + self.distributions.uniform(-1, 1) * self.resamp
        return self._phi


class TwoTinoAsym(Resonant):
    """
    Parameterization for Asymmetric 2:1 resonance

    for the Asymmetric resonance distribute e as a truncated gaussian and select the centre of the
    libration island using the eccentricity.
    """

    def __init__(self, e_mu=0.275, e_sigma=0.06, e_min=0.1, e_max=0.4, j=2, k=1, leading=True, **kwargs):
        super().__init__(j=j, k=k, **kwargs)
        self._phi0 = None
        self.leading = leading
        self.e_mu = e_mu
        self.e_sigma = e_sigma
        self.e_min = e_min
        self.e_max = e_max

    @property
    def e(self):
        """
        Orbital eccentricity picked uniformly between e_min and e_max where e_max must provide q > 22 au.
        """
        if self._e is None:
            # since we are resetting the 'e' distribution make sure to also reset the phi0 distribution
            # here we use phi0 instead of res_centre which is used in the generic resonant class and is constant.
            self._phi0 = None
            # e_max set by q limit.
            e_max = numpy.min(1. - ((22. * units.au) / self.a)).value
            self._e = self.distributions.truncated_normal(self.e_mu,
                                                          self.e_sigma,
                                                          self.e_min,
                                                          e_max)
        return self._e

    @property
    def phi0(self):
        """
        Within the resonance we compute the centre of libration in an eccentricity dependent way.
            - for low-e objects phi is uniform between 140 and 150
            - for high-e we use a hamiltonian expansion (?) to get phi
        """
        if self._phi0 is None:
            # set the phi0 using the polynomial fit provided by Kat Volk.
            self._phi0 = (133.157 -
                          315.272 * self.e +
                          377.082 * self.e ** 2 -
                          0.491667 / self.e) * units.deg
            self._phi0 += (1 - numpy.sin(self.distributions.uniform(0, numpy.pi / 2))) * 20.0 * units.deg
            high_e_cond = self.e > 0.05
            # for low-e orbits a better match is just a constant value.
            self._phi0[~high_e_cond] = self.distributions.uniform(140, 150)[~high_e_cond] * units.deg

            if not self.leading:
                # flip around for trailing island.
                self._phi0 = 360. * units.deg - self._phi0
        return self._phi0

    @property
    def phi(self):
        """
        Phi depends on e,  From Volk/Chan model
        phi52 = (libc + 2*resamp*(ran3(seed) - 0.5))
        """
        if self._phi is None:
            self._phi = self.phi0 + 2*self.resamp*self.distributions.uniform(-0.5, 0.5)
        return self._phi

    @property
    def resamp(self):
        """
        Amplitude of libration around the centre of the libration island.

        Value is picked from a polynomial fit to the resamp vs eccentricity with resamp larger at smaller-e.
        """
        if self._resamp is None:
            # first make resamp appropriate for low-e orbits.
            amp_max = (-403.632 + 9.09917 * self.phi0.to('deg').value - 0.0442498 *
                       self.phi0.to('deg').value ** 2 - 0.0883975 / self.phi0.to('deg').value) * units.deg
            amp_max[self.e < 0.05] = 15 * units.deg
            amp_min = (79.031 * numpy.exp(-(self.phi0.to('deg').value - 121.3435) ** 2 / (2 * 15.51349 ** 2))) * units.deg
            amp_min[self.e < 0.05] = 0 * units.deg
            self._resamp = amp_max - self.distributions.linear(0.25, 1) * (amp_max - amp_min)
            self._resamp[self.e < 0.05] = 15 * units.deg
        return self._resamp


def run(model_filename, detect_filename, characterization_directory, seed, ntrack):
    """
    Using the Resonant defined in ossssim.parametric build model objects and pass them through the survey simulator, save detections

    Args:
        model_filename (str): Name of file to write the generated model to.
        detect_filename (str): Name of file to store detected targets.
        characterization_directory (str): Relative or absolute path to directory on disk where the characterization files are organized
        seed (int): random number seed, specifying allows reproducibility
        ntrack (int): < 0 continue for ntrack iterations;
                      > 0 continue until n_tracked tracked detections;
                      = 0 continue until input exhausted
    """
    ssim = ossssim.OSSSSim(characterization_directory=characterization_directory)

    sym_model = TwoTinoSymmetric(component='Symmetric')
    leading_model = TwoTinoAsym(component='Asym_Lead')
    trailing_model = TwoTinoAsym(leading=False, component='Asym_Trail')

    model_file = ossssim.ModelOutputFile(model_filename)
    model_file.epoch = sym_model.epoch
    model_file.longitude_neptune = sym_model.longitude_neptune
    model_file.colors = ossssim.definitions.COLORS.values()
    model_file.write_header(seed)

    detect_file = ossssim.DetectFile(detect_filename)
    detect_file.epoch = sym_model.epoch
    detect_file.longitude_neptune = sym_model.longitude_neptune
    detect_file.colors = ossssim.definitions.COLORS.values()
    detect_file.write_header(seed)

    n_iter = n_track = n_hits = 0
    while True:
        # pick which model to pass, based on symmetric,trailing/leading ratios
        selection = numpy.random.rand()
        if selection < 0.4:
            model = sym_model
        elif selection < 0.7:
            model = leading_model
        else:
            model = trailing_model
        row = next(model)
        n_iter += 1
        result = ssim.simulate(row, seed=model.seed)
        model_file.write_row(result)
        if result['flag'] > 0:
            n_hits += 1
        if result['flag'] == 4:
            n_track += 1
            detect_file.write_row(result)
        if (0 < ntrack < n_track) or (0 < -ntrack < n_iter):
            break

    detect_file.write_footer(n_iter=n_iter, n_hits=n_hits, n_track=n_track)
    model_file.write_footer(n_iter=n_iter, n_hits=n_hits, n_track=n_track)


def face_down_plot(model_file: str, detect_file: str,
                   detection_table: Table = None) -> None:
    """_
    Plot the detected objects in a face-down plot
    Args:
        model_file: name of file with the base model.
        detect_file: name of file with the detected sources
        detection_table: an astropy Table with helioentric x, y columns
    """
    detect_model = ossssim.ModelFile(detect_file)
    plot = ossssim.RosePlot(detect_model.epoch)
    plot.add_model(detect_model, ms=10, mc='g')
    # noinspection PyArgumentEqualDefault
    drawn_model = ossssim.ModelFile(model_file)
    plot.add_model(drawn_model, mc='k', ms=1, alpha=0.5, sample_size=5000)
    plot.add_pointings('Surveys/OSSOS/pointings.list')
    if detection_table is not None:
        plot.add_detections(detection_table)
    plot.add_galactic_plane()
    plot.show()


if __name__ == '__main__':

    # Goal is to model the OSSOS 2:1 detections.  First we load all the OSSOS detections from the
    # summary file and determine how many 2:1 objects there were.

    data = Table.read(os.path.join(SURVEYS_DIR,'OSSOSv11/ObsSummary/OSSOS/OSSOS.CDS')
                      , format='cds')
    cond = (data['cl'] == 'res') & (data['j'] == 2) & (data['k'] == 1)

    # now run a survey simulation.
    run('TwotinoModel.dat',
        'TwotinoDetect.dat',
        os.path.join(SURVEYS_DIR,'OSSOSv11/ObsSummary/OSSOS'),
        123456789,
        55)

    # confirm this looks like the OSSOS detections using rose plot.
    face_down_plot('TwotinoModel.dat',
                   'TwotinoDetect.dat',
                   detection_table=data)
