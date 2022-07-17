"""
Abstract class for building parametric models from.
"""
from abc import ABC, abstractmethod

import numpy
from astropy import units
from astropy.table import QTable
from astropy.units import Quantity

from . import distributions
from .core import OSSSSim


class HDistribuion(ABC):
    """
    Provide a class to describe the size distribiution of the object being simulated.
    """
    def __init__(self, func, **kwargs):
        """
        Args:
            func (func): the function that will be used for the H-distribution.
            kwargs (dict): A dictionary of parameters used by the H distribution function.
        """
        self.params = kwargs
        self.func = func

    def __call__(self):
        return self.func(**self.params)


class NonResonant(ABC):
    """
    This abstract class defines methods needed to build a parametric Outer Solar System model for use as a model input for OSSSSIM
    """

    def __init__(self,
                 size=10**6,
                 seed=123456789,
                 epoch=2456293.5 * units.day):
        """
        Setup the boundaries of the simulation.  size and seed are used to initialize a dist_utils.Distribution class.

        Args:
            size (int): Determines size of the arrays to be generated (default=10^6).
            seed (int): Initialize distributions with this seed to enable reproducible models.

        """
        if seed is None:
            seed = numpy.random.randint(1, 999999999)
        self.seed = seed
        self.size = size
        self._epoch = epoch
        self.distributions = distributions.Distributions(self.seed, self.size)
        self.h_distribution = HDistribuion(self.distributions.power_knee_divot,
                                           **dict([('alpha_b', 1.1),
                                                   ('alpha_f', 0.4),
                                                   ('h_break', 7.5),
                                                   ('h_max', 11),
                                                   ('h_min', 1)]))

        # these are protected variables
        self._iter = None  # holds an iterator basd on the iterable QTable that holds element values
        self._targets = None  # holds the value of the property target
        self._comp = "Cls"

    @property
    @abstractmethod
    def a(self):
        """
        J2000 Heliocentric semi-major axis.
        """
        pass

    @property
    @abstractmethod
    def e(self):
        """
        J2000 Heliocentric eccentricity
        """
        pass

    @property
    @abstractmethod
    def inc(self):
        """
        J2000 heliocentric inclination of orbit
        """
        pass

    @property
    def node(self):
        """
        Return uniformly distributed nodes.
        """
        return self.distributions.uniform(0, 2*numpy.pi) * units.rad

    @property
    def peri(self):
        """
        Distribute peri centre to obey the phi/M/_longitude_neptune constraints.
        """
        return self.distributions.uniform(0, 2*numpy.pi) * units.rad

    @property
    def M(self):
        """
        Return uniformly distributed mean anomalies
        """
        return self.distributions.uniform(0, 2*numpy.pi) * units.rad

    @property
    def H(self):
        """A distribution of H values"""
        # define the default size distribution parameters.
        return self.h_distribution()

    @property
    def epoch(self):
        """
        Epcoh, in Julian Days, of the orbital elements.
        """
        if self._epoch is not None:
            self._epoch = self.distributions.constant(self._epoch) * units.day
        return self._epoch

    @property
    def resamp(self):
        """
        Retuns resamp a 0.0 as this is not a resonant orbit.
        """
        return self.distributions.uniform(0.0)

    @property
    def j(self):
        """
        Return j of (j/k) as 0 as this is not a resonant orbit
        """
        return self.distributions.uniform(0)

    @property
    def k(self):
        """
        Return k of (j:k) as 0 as this is not a resonant orbit
        """
        return self.distributions.uniform(0)

    @property
    def comp(self):
        """
        Label for the component being generated.
        """
        return self._comp * self.size

    @comp.setter
    def comp(self, value):
        self._comp = value

    @property
    def lc_gb(self):
        """
        Opposition surge effect as define in Bowell
        """
        return self.distributions.constant(OSSSSim.LIGHT_CURVE_PARAMS['gb']) * units.mag

    @property
    def lc_phase(self):
        """
        Phase of lightcurve at self.epoch
        """
        return self.distributions.constant(OSSSSim.LIGHT_CURVE_PARAMS['phase']) * units.rad

    @property
    def lc_period(self):
        """
        period of lightcurve
        """
        return self.distributions.constant(OSSSSim.LIGHT_CURVE_PARAMS['period']) * units.day

    @property
    def lc_amplitude(self):
        """
        peak-to-peak amplitude of lightcurve
        """
        return self.distributions.constant(OSSSSim.LIGHT_CURVE_PARAMS['amplitude']) * units.mag

    @property
    def colors(self):
        """
        colors of objects expressed as a list of list:

        [ [ (g-x), (g-r), (i-x), (z-x), (u-x), (V-x), (B-x), (R-x), (I-x), (w-x) ], ... ]
        """

        g_x = self.distributions.normal(0.7, 0.2) * units.mag  # mu 0.7 and std 0.2 (from CFEPS data)
        r_x = self.distributions.constant(0) * units.mag
        i_x = self.distributions.constant(0) * units.mag
        z_x = self.distributions.constant(0) * units.mag
        u_x = self.distributions.constant(0) * units.mag
        V_x = self.distributions.constant(0) * units.mag
        B_x = self.distributions.constant(0) * units.mag
        R_x = self.distributions.constant(-0.1) * units.mag
        I_x = self.distributions.constant(0) * units.mag
        w_x = self.distributions.constant(0) * units.mag

        return numpy.column_stack((g_x, r_x, i_x, z_x, u_x, V_x, B_x, R_x, I_x, w_x))

    def _generate_targets(self):
        """
        Generate the orbit elements and properties of a list of targets to be passed to the simulator.

        Expected to return a QTable or dictionary.  If a QTable then len of table should be self.size.
        If dictionary then each dictionary key should point to a list of length self.size.

        All values should stored as Quantity objects to allow conversion to desired units before passing to the SurveySubsF95.detos1

        Must define at least {'a': [], 'e': [], 'inc': [], 'node': [], 'peri': [], 'M': [], 'H': []}  see OSSSSIM.simulate for full
        list of keys that can be returned.

        Returns:
            (QTable or dict): set of Quantity objects describing targets.
        """
        return QTable([self.a, self.e, self.inc, self.node, self.peri, self.M, self.epoch,
                       self.H, self.lc_gb, self.lc_phase, self.lc_period, self.lc_amplitude,
                       self.resamp, self.colors, self.comp, self.j, self.k],
                      names=['a', 'e', 'inc', 'node', 'peri', 'M', 'epoch',
                             'H', 'lc_gb', 'lc_phase', 'lc_period', 'lc_amplitude',
                             'resamp', 'colors', 'comp', 'j', 'k'])

    @property
    def targets(self):
        """
        A table of length self.size holding the generated distribution of targets.

        """
        if self._targets is None:
            self._targets = self._generate_targets()
        return self._targets

    @targets.setter
    def targets(self, value):
        """
        set targets to Table stored in value.

        Args:
            value (Table or None): value to set targets to.
        """
        if not isinstance(value, QTable) or value is not None:
            raise ValueError(f"Attempted to set targets table to something that isn't a table: {value}")
        self._targets = value

    @property
    def iter(self):
        """
        An iterator on self.targets table
        """
        if self._iter is None:
            self._iter = iter(self.targets)
        return self._iter

    def __iter__(self):
        return self

    def __next__(self):
        """
        return the next row from the orbits table via the `iterator` on the table.

        If we hit the end of orbits table then call draw_distribution to refresh the table.
        """
        try:
            row = next(self.iter)
        except StopIteration:
            # Clear the targets table so a new distribution will be generated.
            self.targets = None
            row = next(self.iter)

        return row


class Resonant(NonResonant):
    """
    This class defines methods needed to build a parametric Outer Solar System model for use as a model input for OSSSSIM
    """

    def __init__(self,
                 size=10**6,
                 seed=123456789,
                 component='Res',
                 longitude_neptune=333.5078 * units.deg,
                 epoch_neptune=2456293.5 * units.day,
                 j=2,
                 k=2,
                 res_amp_low=0*units.deg,
                 res_amp_mid=5*units.deg,
                 res_amp_high=10*units.deg,
                 res_centre=0*units.deg):
        """
        Setup the boundaries of the simulation.  size and seed are used to initialize a dist_utils.Distribution class.

        _longitude_neptune and epoch_neptune are stored in self._longitude_neptune and self.epoch_neptune for use in
        self._generate_targets.

        See examples/paramertic.py for an example implementation.

        Args:
            size (int): Determines size of the arrays to be generated (default=10^6).
            seed (int): Initialize distributions with this seed to enable reproducible models.
            longitude_neptune (Quantity): The ecliptic longitude of Neptune at epoch_neptune of elements.
            epoch_neptune (Quantity): The epoch_neptune of the given longitude of Neptune.

        """
        super(Resonant).__init__(size=size, seed=seed, epoch=epoch_neptune)

        self.longitude_neptune = longitude_neptune  # Neptune's mean longitude on 1 Jan 2013
        self.a_neptune = 30.07 * units.au
        self.epoch_neptune = epoch_neptune  # this epoch_neptune corresponds to _longitude_neptune
        self._j = j
        self._k = k
        self.res_amp_low = res_amp_low
        self.res_amp_high = res_amp_high
        self.res_amp_mid = res_amp_mid
        self.res_centre = res_centre
        self._comp = component.replace(" ","_")

    @property
    def a(self):
        """
        J2000 Heliocentric semi-major axis sampled as +/- 0.5 from the resonance semi-major axis value.
        """
        a0 = (self.a_neptune ** (3 / 2) * self.j / self.k) ** (2 / 3)
        amin = a0 - 0.5 * units.au
        amax = a0 + 0.5 * units.au
        return self.distributions.uniform(amin, amax)

    @property
    def e(self):
        """
        Set the maximum value of 'e' based on the peri-center location of Neptune, minimum value set to 0.02 then randomly sample this
        range of e.
        """
        e_max = (1 - self.a_neptune / self.a)
        e_min = 0.02
        return self.distributions.uniform(e_min, e_max)

    @property
    def inc(self):
        """
        Distribute the inclinations based on Brown 2001 functional form.
        """
        return self.distributions.truncated_sin_normal(0*units.degree.to('rad'),
                                                       11*units.degree.to('rad'),
                                                       0*units.degree.to('rad'),
                                                       40*units.degree.to('rad') * units.rad)

    @property
    def peri(self):
        """
        Distribute peri centre to obey the phi/M/_longitude_neptune constraints.
        """
        return (1 / 1 * (self.phi - 2 * self.M) - self.node + self.longitude_neptune) % (2 * numpy.pi)

    @property
    def phi(self):
        """
        Compute the phi, librartion centre from the resonance centre and sampling the resonance amplitude via sin() weighting.
        """
        return self.res_centre + numpy.sin(self.distributions.uniform(0, 2*numpy.pi))*self.resamp

    @property
    def resamp(self):
        """
        amplitude of the distribution of resonance centres around libration centre, used in self.phi
        """
        return self.distributions.linear_peak(self.res_amp_low,
                                              self.res_amp_mid,
                                              self.res_amp_high)

    @property
    def j(self):
        """
        The j component of the j:k resonance definition
        """
        return self.distributions.constant(self._j)

    @j.setter
    def j(self, value):
        """
        Set the j of j/k resoance discription.
        """
        self._j = value

    @property
    def k(self):
        """
        The k component of the j:k resoance definition
        """
        return self.distributions.constant(self._k)

    @k.setter
    def k(self, value):
        """
        Set the k of j/k resoance discription.
        """
        self._k = value
