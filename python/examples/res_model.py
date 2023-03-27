"""
An example of a 2:1 model for use with OSSSSIM.simulate.  This is the OSSOS++ 2:1 model from Chen and Volk.

To work in the python mode here we define two new classes derived from the Resonance.  The Resonance base class handles computing
peri from the libration centre and by default will distributed the libration centre using a 'low, mid, high' set of values that define
a triangle distribution (starts at 0 at low, rises to 1 at mid, back to 0 at high) for the amplitude distribution.  That was CFEPS.

For OSSOS we have a more complex set computations where the libration centre (phi0) and the libration amplitude (resamp) about that
 centre depend on the eccentricity of orbit and type (Symmetric or Asymmetric) of the resonance.
"""
import os
import sys

import numpy
from astropy import units
from astropy.table import Table
from matplotlib import pyplot as plt

import ossssim
from funcs import variably_tapered2
from ossssim.models import Resonant
from ossssim.plotter import RosePlot

# From JPL Horizons
# Ecliptic BaryCentre location of New Horizons
nhtraj = {'jd': 2460563.500000000,
          'x': 1.853275545830267E+01 * units.au,
          'y': -5.684640535274954E+01 * units.au,
          'z': 2.076736425369225E+00 * units.au,
          'vx': 3.089351377914578E-03 * units.au / units.day,
          'vy': - 7.264635908893537E-03 * units.au / units.day,
          'vz': 2.856560416694524E-04 * units.au / units.day}

nhtraj['r'] = (nhtraj['x']**2+nhtraj['y']**2 + nhtraj['z']**2)**0.5
nhtraj['vel'] = (nhtraj['vx']**2+nhtraj['vy']**2 + nhtraj['vz']**2)**0.5

print(f"NH: distance: {nhtraj['r'].to('au')} velocity: {nhtraj['vel'].to('km/s')}")

def H_cfd(H):
    alpha_SI = 5*0.4/3
    beta_SI = 0.42
    Ho = -2.6
    Hb = 8.1
    return variably_tapered2(H, Ho, Hb, alpha_SI, beta_SI)


class Symmetric(Resonant):
    """
    Build a Symmetric class from the resonant base class, j/k to 2/1.
    """
    def __init__(self, q_c=0.07, q_w=0.30, j=2, k=1, sigma_i=14.5, H_min=4.5, H_max=17.0, **kwargs):
        """
        Symmetric, libration centers are all 180
        :param q_c: center of q distribution
        :param q_w: width of q distribution
        :param j: neptune orbits
        :param k: tno orbits
        
        pick libration amplitudes uniformly from 80-160
        pick e uniformly from 0.05 -> 0.3 (Chiang&Jordan02)

                res_amp_low=20*units.degree, 
                res_amp_mid=95*units.degree,
                res_amp_high=130*units.degree,

        """
        # amplitude values are from Gladman et al. 2016
        super().__init__(j=j, k=k, 
                **kwargs)
        self._q = None
        self.q_c = q_c
        self.q_w = q_w
        self.sigma_i = sigma_i
        self.res_centre = 180.0 * units.deg
        self.H_min = H_min
        self.H_max = H_max


    @property
    def inc(self):
        """
        Distribute the inclinations based on Brown 2001 functional form.
        """
        if self._inc is None:
            self._inc = self.distributions.truncated_sin_normal(0,
                                                                numpy.deg2rad(self.sigma_i),
                                                                0,
                                                                numpy.deg2rad(45)) * units.rad
        return self._inc

    @property
    def e(self):
        """
        Set values of 'e' by sampling 'q' and the setting 'e'
        """
        if self._e is None:
            # self._e = self.distributions.uniform(0.3,0.33)
            # return self._e
            # max is set by q but also limited by users choice of e_max.
            res_a = 29.9*((self.j[0]/self.k[0])**(2/3))
            q = self.distributions.truncated_normal(self.q_c, self.q_w, res_a*(1-0.8), res_a*(1-0.001))
            self._e = 1 - q/res_a
        return self._e

    @property
    def off_resamp(self):
        """
        Uniformly distributed across the width
        """
        if self._resamp is None:
            # self._resamp = self.distributions.uniform(0., 1.) * units.deg
            # return self._resamp
            self._resamp = self.distributions.uniform(80., 160.) * units.deg
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

    @property
    def H(self):
        """
        Generate an H distribution following Kavelaars et al. (2021).  
        """
        if self._H is None:
            fp = numpy.arange(self.H_min, self.H_max, 0.1)
            xp = H_cfd(fp)
            xp = xp/xp[-1]
            x = self.distributions.rnd_gen.uniform(xp[0], xp[-1], self.distributions.size)
            self._H = numpy.interp(x, xp, fp)
        return self._H


class Asymmetric(Symmetric):
    """
    Parameterization for Asymmetric n:1 resonance  which is taken from Volk

    for the Asymmetric resonance distribute e as a truncated gaussian and select the centre of the
    libration island using the eccentricity.
    """

    def __init__(self, leading=True, **kwargs):
        super().__init__(**kwargs)
        self.leading = leading
        self._phi0 = None

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
        Phi0 depends on e,  From Volk/Chan model
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
            # self._resamp = self.distributions.uniform(0., 1.) * units.deg
            # return self._resamp
            # first make resamp appropriate for low-e orbits.
            amp_max = (-403.632 + 9.09917 * self.phi0.to('deg').value - 0.0442498 *
                       self.phi0.to('deg').value ** 2 - 0.0883975 / self.phi0.to('deg').value) * units.deg
            amp_max[self.e < 0.05] = 15 * units.deg
            amp_min = (79.031 * numpy.exp(-(self.phi0.to('deg').value - 121.3435) ** 2 / (2 * 15.51349 ** 2))) * units.deg
            amp_min[self.e < 0.05] = 0 * units.deg
            self._resamp = amp_max - self.distributions.linear(0.25, 1) * (amp_max - amp_min)
            self._resamp[self.e < 0.05] = 15 * units.deg
        return self._resamp


def LongLat(x, y, z):
    """
    Compute Spherical long/lat of cartesian position.
    """
    r = (x**2 + y**2 + z**2)**0.5
    long = numpy.arctan2(y, x)
    lat = numpy.arcsin(z / r)
    return long, lat, r


def run(model_filename, params, seed, H_max=14.0):
    """
    Using the Resonant defined in ossssim.parametric build model objects and pass them through the survey simulator, save detections

    Args:
        model_filename (str): Name of file to write the generated model to.
        params (str): name of file with the resonance parameters
        seed (int): an integer value to seed the random generator
        H_max(float): the maximum H value to simulate.

    params file should have columns 'j', 'k', 'q_c', 'q_w', 'sigma_i', 'MedianPop'

    MedianPop is expected to be to magnitude H=8.66 (OSSOS Standard value)
    P(i) ~ sin(i) * exp(-i**2/2*sigma_i**2)
    P(q) ~ exp(-(q-q_c)**2/2*q_w**2)
    """

    model_file = ossssim.ModelOutputFile(model_filename)
    model_file.colnames = ['a', 'e', 'inc', 'node', 'peri', 'M', 'H', 'x', 'y', 'z', 'delta',
                           'vx', 'vy', 'vz', 'delta_v', 'dt', 'comp']
    # model_file.colors = ossssim.definitions.COLORS.values()
    model_file.longitude_neptune=6.876 * units.rad
    model_file.epoch=2456839.5 * units.day
    model_file.write_header(seed)

    res_table = Table.read(params, format='ascii')
    H_model = 8.66
    N1 = H_cfd(H_model)
    N2 = H_cfd(H_max)
    print(f"Scaling up by {N2/N1}")
    for row in res_table:
        models = {}
        j=int(row['j'])
        k=int(row['k'])
        q_c = float(row['q_c'])
        q_w = float(row['q_w'])
        sigma_i = float(row['sigma_i'])
        comp=f"Res-{j}:{k}"
        N = int(row['MedianPop']*N2/N1)
        print(f"Generating {N} sources for row:\n{row}")
        components = ['sym', 'leading', 'trailing']
        selections = {components[0]: 1.0,
                      components[1]: 2.0,
                      components[2]: 2.0}
        models[components[0]] = Symmetric(j=j, k=k, component=comp,
                   q_c=row['q_c'], q_w=row['q_w'], sigma_i=sigma_i, H_max=H_max, res_amp_low=20.0*units.degree,
                   res_amp_high=160*units.degree, res_amp_mid=95*units.degree, size=min(N, 1000000))
        if k == 1:
            selections = {components[0]: 0.3,
                          components[1]: 0.3 + 0.35,
                          components[2]: 1}
            models[components[1]] = Asymmetric(j=j, k=k, component=comp,
                                               q_c=q_c, q_w=q_w, sigma_i=sigma_i, H_max=H_max, size=min(N, 1000000))
            models[components[2]] = Asymmetric(leading=False, j=j, k=k, component=comp,
                                               q_c=q_c, q_w=q_w, sigma_i=sigma_i, H_max=H_max, size=min(N, 1000000))

        niter = 0
        while niter<N:
            niter += 1
            selection = numpy.random.rand()
            # print(components, selection, selections)
            component = list(filter(lambda x: selection<selections[x], components))[0]
            particle = next(models[component])
            delta = (particle['x']**2 + particle['y']**2 + particle['z']**2)**0.5
            particle_v = (particle['vx']**2 + particle['vy']**2 + particle['vz']**2)**0.5
            if delta < nhtraj['r']:
                continue
            dpos = [particle[c] - nhtraj[c] for c in ['x', 'y', 'z']]
            dvel = [nhtraj[c] - particle[c] for c in ['vx', 'vy', 'vz']]
            long, lat, r = LongLat(dpos[0].to('au'), dpos[1].to('au'), dpos[2].to('au'))
            vlong, vlat, vel = LongLat(dvel[0].to('au/year'), dvel[1].to('au/year'), dvel[2].to('au/year'))
            dt = r / vel
            sep = r*numpy.arccos(numpy.sin(vlat)*numpy.sin(lat)+numpy.cos(vlat)*numpy.cos(lat)*numpy.cos(vlong-long))/units.rad
            delta_v = numpy.fabs(sep / dt)
            if delta_v < 3.00*units.km/units.second:
                tdict = {}
                for col in particle.colnames:
                    tdict[col] = particle[col]
                tdict['delta_v'] = delta_v
                tdict['dt'] = dt
                tdict['delta'] = delta.to('au').value
                model_file.write_row(tdict)

    model_file.write_footer(n_iter=0, n_hits=0, n_track=0)


def face_down_plot(model_file: str) -> None:
    """_
    Plot the detected objects in a face-down plot
    Args:
        model_file: name of file with the base model.
    """
    drawn_model = ossssim.ModelFile(model_file)
    plot = RosePlot(drawn_model.epoch)
    plot.add_model(drawn_model, mc='k', ms=1, alpha=0.5, sample_size=5000)
    plt.savefig(f'{os.path.splitext(model_file)[0]}.pdf')


def main():
    """
    Execute this module as a script.
    """
    # Goal is to model the OSSOS resonance detections given a file with parameters for those resonances.
    # e.g. from Crompvoets et al. (2021)

    # now run a survey simulation.
    params = sys.argv[1]
    H_max = float(sys.argv[2])
    outfile=f"{os.path.splitext(params)[0]}_Model.dat"
    print(f"Saving results to {outfile}")
    if not os.access(outfile, os.R_OK):
        run(outfile, params, 123456789, H_max=H_max)

    # confirm this looks like the OSSOS detections using rose plot.
    face_down_plot(outfile)


if __name__ == '__main__':
    main()


