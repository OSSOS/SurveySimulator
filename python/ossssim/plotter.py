"""
Some methods to aid making orbit plots from simulator outputs
"""
import logging

import numpy
from astropy import units
from astropy.coordinates import SkyCoord
from astropy.time import Time
from astropy.units import Quantity
from astroquery import jplhorizons
from matplotlib import pyplot as plt, cycler
from matplotlib import rcParams
from matplotlib.ticker import MultipleLocator
from numpy.random import default_rng

from . import definitions
from .models import ModelFile
from .models import Parametric

np = numpy

# setup the plotting Look and Feel.
rcParams['font.size'] = 12  # good for posters/slides
rcParams['patch.facecolor'] = (0.4, 0.7607843137254902, 0.6470588235294118)  # brewer2mpl 'Set2' 'qualitative' colors
rcParams['figure.figsize'] = (10, 10)
rcParams['figure.dpi'] = 150
# this is the puor color cycle from brewer2mpl
rcParams['axes.prop_cycle'] = cycler('color', [(0.4980392156862745, 0.23137254901960785, 0.03137254901960784),
                                               (0.7019607843137254, 0.34509803921568627, 0.023529411764705882),
                                               (0.8784313725490196, 0.5098039215686274, 0.0784313725490196),
                                               (0.9921568627450981, 0.7215686274509804, 0.38823529411764707),
                                               (0.996078431372549, 0.8784313725490196, 0.7137254901960784),
                                               (0.9686274509803922, 0.9686274509803922, 0.9686274509803922),
                                               (0.8470588235294118, 0.8549019607843137, 0.9215686274509803),
                                               (0.6980392156862745, 0.6705882352941176, 0.8235294117647058),
                                               (0.5019607843137255, 0.45098039215686275, 0.6745098039215687),
                                               (0.32941176470588235, 0.15294117647058825, 0.5333333333333333),
                                               (0.17647058823529413, 0.0, 0.29411764705882354)])
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Sofia Pro', 'Tahoma']
ALMOST_BLACK = '#262626'


class TimeSeriesPlot:
    """
    Plot the elements of model as function of index in the list, aids in diagnosis the distribution of elements.
    """

    def __init__(self, model: (ModelFile or Parametric)) -> None:
        """Given an input model setup methods to make plots.
        Args:
            model: a model object from the ossssim package."""
        self.model = model
        self.fig = plt.figure(figsize=(8, 15))

    def plot(self, variables: list = None) -> None:
        """
        Plot model elements as the value vs position in the array (index).

        Args:
            variables ([] or None): a list of columns to plot. e.g. (['a', 'e', 'inc']) default:None plots all columns in model.

        To see what columns are available look at self.model.colnames

        Intended for diagnostic purposes.
        """

        if variables is None:
            variables = self.model.targets.colnames
        nx = len(variables)//2
        ny = 2
        for i, column in enumerate(variables):
            ax = self.fig.add_subplot(nx, ny, i)
            if column not in self.model.targets.colnames:
                logging.warning(f"Could not plot {column} as does not appear input model.")
                pass

            values = self.model.targets[column]
            if isinstance(values[0], (list, numpy.ndarray)):
                for idx in range(len(values[0])):
                    ax.plot(values[:, idx].to(definitions.colunits[column]).value,
                            color='k', marker='o', linestyle='none', linewidth=2, markersize=1)
            else:
                ax.plot(values.to(definitions.colunits[column]).value,
                        color='k', marker='o', linestyle='none', linewidth=2, markersize=1)
            ax.ylabel(f"{column} ({definitions.colunits[column]})")
        plt.show()


class RosePlot:
    """
    Plot objects of model as face-down view of the outer solar system.
    Construct a plot of solar orbits (based on OSSSSim ModelFile and Characterization and OSSOS
    formatted detection lists) in a top-down view.

    We have called this the 'Rose Plot' as the blocks drawn on the top-down view give a sort of 'rose petal' look.. it's a stretch.

    """

    def __init__(self, epoch: float, outer_edge=85 * units.au, inner_edge=10 * units.au) -> None:
        """
        Plot the TNO discoveries on a top-down Solar System showing the position of Neptune and model TNOs.

        Discoveries are plotted at their time of discovery according to the positions in detection file, formatted like OSSOS.CDS file

        Coordinates are polar RA, radial axis is in AU.

        This is not a 'projection' of the particles into any common plane.

        Each wedge is a different latitude above the plane, but inside each wedge it is heliocentric distance vs RA.
        If someone asks we know this is NOT (for example) a projection of each object down into the ecliptic.

        TODO: Make galactic wedge more realistic.
        """
        self.epoch = epoch
        self._longitude_neptune = None
        self.frame = 'heliocentrictrueecliptic'
        self.outer_edge = outer_edge.to('au').value
        self.inner_edge = inner_edge.to('au').value
        self.fig = plt.figure(figsize=(8, 8))
        rect = [0.0725, 0.0725, 0.85, 0.85]  # the plot occupies not all the figure space
        self.ax1 = self.fig.add_axes(rect, polar=True, frameon=False)  # theta (RA) is zero at E, increases anticlockwise
        self.ax1.set_aspect('equal')

        self.ax1.set_rlim(0, self.outer_edge)
        rings = range(0, 100, 15)
        ring_labels = ["", ""].extend([f"{x:3d} au" for x in rings])
        self.ax1.set_rgrids(rings, ring_labels, angle=100, alpha=0.45)
        self.ax1.yaxis.set_major_locator(MultipleLocator(25))
        self.ax1.xaxis.set_major_locator(MultipleLocator(numpy.deg2rad(15)))  # every 2 hours
        self.ax1.grid(axis='x', color='k', linestyle='--', alpha=0.2)
        x_tick_labels = []
        lon = numpy.arange(0, 360, 30) * units.deg
        lat = numpy.zeros(len(lon)) * units.deg
        dist = numpy.ones(len(lon)) * 45 * units.au
        coord = SkyCoord(lon, lat, distance=dist, obstime='2000-01-01', frame='heliocentrictrueecliptic').transform_to('icrs')

        for label_values in coord.ra.hour:
            lv = int(numpy.round(label_values))
            x_tick_labels.append('')
            x_tick_labels.append(f"{lv}h")
        self.ax1.set_xticklabels(x_tick_labels, color='b', alpha=0.6)

    @property
    def longitude_neptune(self):
        """
        Return the longitude of Neptune based on the current epoch
        """
        if self._longitude_neptune is None:
            neptune = jplhorizons.Horizons(899, epochs=self.epoch, location='568')
            self._longitude_neptune = neptune['EcLon']
        return self._longitude_neptune

    def add_pointings(self, pointing_filename, color='b', alpha=0.1, label=False):
        """
        Read in an OSSOS pointing file and add the blocks in that file to the RosePlot
        """
        names = []
        with open(pointing_filename, 'r') as file_obj:
            while True:
                line = file_obj.readline()
                if len(line) == 0:
                    break
                if line.startswith('#'):
                    continue
                attributes = line.split()
                if label:
                    name = attributes[-1][0:2]
                    if name in names:
                        name = None
                    names.append(name)
                else:
                    name = None
                pos = SkyCoord(attributes[2],
                               attributes[3],
                               unit=('hour', 'deg', 'au'),
                               obstime='2000-01-01',
                               distance=44).transform_to(self.frame)
                if not (-10 < pos.lat.degree < 10) :
                    color='y'
                if 'poly' in attributes[0]:
                    n_vertices = int(attributes[1])
                    min_ra = max_ra = None
                    for idx in range(n_vertices):
                        vertices = [float(x) for x in file_obj.readline().split()]
                        if min_ra is None or min_ra > vertices[0]:
                            min_ra = vertices[0]
                        if max_ra is None or max_ra < vertices[0]:
                            max_ra = vertices[0]
                    width = (max_ra - min_ra) * units.deg
                else:
                    width = float(attributes[0]) * units.deg
                self.add_block(pos.lon, width, color=color, alpha=alpha, label=name)

    def add_block(self, ra_cen: Quantity, width: Quantity, color='b',
                  alpha=0.1, label=None) -> None:
        """
        Add a 'wedge' to the polar plot to show the location of a block.
        Args:
            ra_cen: the Right Ascension of the centre of the block
            width: the full width of the block
            color: face color to fill in the block
            alpha: transparency of block
            label: name of the block
        """
        self.ax1.bar(ra_cen.to('rad').value,
                     self.outer_edge,
                     linewidth=0.1,
                     width=width.to('rad').value,
                     bottom=self.inner_edge,
                     zorder=0,
                     color=color,
                     alpha=alpha)
        if label is not None:
            self.ax1.annotate(label, (ra_cen.to('rad').value,
                              self.outer_edge + 3.),
                              size=25, color=ALMOST_BLACK)

    def add_galactic_plane(self, minimum_latitude=22 * units.deg):
        """
        Put a grey bar over the area that has the galactic plane contamination,
        ie. the ecliptic plane contains galactic latitude < minimum_latitude
        """
        lon = numpy.arange(0, 360, 0.1)*units.deg
        lat = 0*lon
        coords = SkyCoord(lon,
                          lat,
                          distance=40 * units.au,
                          frame=self.frame,
                          obstime=Time('2000-01-01')).transform_to('galactic')
        start_arc = end_arc = None
        for coord in coords:
            if minimum_latitude > coord.b > -minimum_latitude:
                if start_arc is None:
                    start_arc = coord.transform_to(self.frame)
            if (coord.b > minimum_latitude or coord.b < -minimum_latitude) and start_arc is not None:
                end_arc = coord.transform_to(self.frame)
            if end_arc is not None and start_arc is not None:
                end_lon = end_arc.lon.radian
                start_lon = start_arc.lon.radian
                if end_lon < start_lon:
                    end_lon += 2*numpy.pi
                plane = (end_lon + start_lon)/2.0
                width = end_lon - start_lon
                self.ax1.bar(plane, self.outer_edge, width=width, color=ALMOST_BLACK, linewidth=0, alpha=0.2)
                self.ax1.annotate('galactic plane', (plane, self.outer_edge - 15), size=10, color='k', alpha=0.45)
                end_arc = start_arc = None

    # No optional on these just yet
    def add_planets(self):
        """
        Use Horizons to lookup the positions of the planets at epoch and then plot those locations.
        """
        ids = {'Jupiter': 599, 'Saturn': 699, 'Uranus': 799,
               'Neptune': 899, 'Pluto': 999}
        fc = ALMOST_BLACK
        for planet_name in ['Jupiter', 'Saturn', 'Uranus', 'Neptune', 'Pluto']:
            planet = jplhorizons.Horizons(ids[planet_name],
                                          location='568',
                                          epochs=self.epoch)
            eph = planet.ephemerides()
            alpha = 0.7
            size = 20
            if planet_name == 'Pluto':
                alpha = 0.35
                size = 10
            self.ax1.scatter(eph['EclLon'].to('rad').value, eph['r'].to('au').value,
                             marker='o',
                             s=size,
                             facecolor=fc,
                             edgecolor=fc,
                             alpha=alpha)

    def add_detections(self, detection_table):
        """
        Taking a OSSOS/CFEPS detection list in CDS format an plot real detections.

        Args:
            detection_table (Table): astropy table to add to face_down plot must have columns 'x' and 'y' in 'au'
        """
        # noinspection SpellCheckingInspection
        coords = SkyCoord(detection_table['RAdeg'],
                          detection_table['DEdeg'],
                          obstime='2000-01-01',
                          distance=detection_table['dist'],
                          frame='icrs').transform_to(self.frame)
        self.ax1.scatter(coords.lon.to('rad').value,
                         coords.distance.to('au').value,
                         s=5,
                         c='c')

    def add_scale_rings(self, radii: list = [10, 30, 50, 100]) -> None:
        """
        Add scale rings to the plot.

        radii: the radii of rings to plot to provide scale to the plot.
        """

        theta = numpy.arange(0, 2*numpy.pi, 2*numpy.pi/1000)
        for guide_circle in radii:
            r = numpy.ones(len(theta))*guide_circle
            self.ax1.plot(theta, r, ls=':')
            self.ax1.annotate(f"{guide_circle} au", (0, guide_circle+1),
                              color='b', horizontalalignment='center',
                              verticalalignment='center_baseline')

    def add_model(self, model: ModelFile, mc: str = 'k', ms: float = 1.,
                  sample_size: int = None, alpha: float = 1.0) -> None:
        """
        Make a face-down plot of the solar system for this model.
        """
        state = model.cartesian
        coord = SkyCoord(state['x'], state['y'], state['z'], representation_type='cartesian',
                         frame='heliocentrictrueecliptic', obstime='2000-01-01').transform_to(self.frame)

        if sample_size is None:
            choice = numpy.arange(len(coord))
        else:
            rng = default_rng()
            choice = rng.integers(0, len(model.targets), sample_size)

        self.ax1.plot(coord.lon.to('rad').value[choice],
                      coord.distance.to('au').value[choice],
                      f'.{mc}',
                      markersize=ms,
                      alpha=alpha)

    @staticmethod
    def show() -> None:
        """
        Display the plot.
        """
        plt.show()
