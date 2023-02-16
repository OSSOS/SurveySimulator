__author__ = 'Michele Bannister'
import math

import ephem
import matplotlib.pyplot as plt
import numpy as np
from astropy import units
from astropy.units import Quantity
from matplotlib import rcParams
from matplotlib.ticker import MultipleLocator

class RosePlot:
    """
    Construct a plot of solar orbits (based on OSSSSim ModelFile and Characterization and OSSOS
    formatted detection lists) in a top-down view.

    We have called this the 'Rose Plot' as the blocks drawn on the top-down view give a sort of 'rose petal' look.. it's a stretch.
    """

    def __int__(self, outer_edge=85*units.au, inner_edge=10*units.au):
        """
        Plot the TNO discoveries on a top-down Solar System showing the position of Neptune and model TNOs.

        Discoveries are plotted at their time of discovery according to the positions in detection file, foramted like OSSOS.CDS file

        Coordinates are polar RA, radial axis is in AU.

        This is not a 'projection' of the particles into any common plane.

        Each wedge is a different latitude above the plane, but inside each wedge it is heliocentric distance vs RA.
        If someone asks we know this is NOT (for example) a projection of each object down into the ecliptic.

        TODO: Make galactic wedge more realistic.
        """


        self.outer_edge = outer_edge
        self.inner_edge = inner_edge
        self.fig = plt.figure(figsize=(6, 6))
        rect = [0.01, 0.01, 0.95, .95]  # the plot occupies not all the figspace
        self.ax1 = self.fig.add_axes(rect, polar=True, frameon=False)  # theta (RA) is zero at E, increases anticlockwise
        self.ax1.set_aspect('equal')

        self.ax1.set_rlim(0, outer_edge)
        rings = range(20, 100, 20)
        ring_labels = ["", ""].extend([ f"{x:3d} au" for x in rings])
        self.ax1.set_rgrids(rings, ring_labels, angle=90, alpha=0.45)
        self.ax1.yaxis.set_major_locator(MultipleLocator(20))
        self.ax1.xaxis.set_major_locator(MultipleLocator(math.radians(15)))  # every 2 hours
        self.ax1.grid(axis='x', color='k', linestyle='--', alpha=0.2)
        xticklabels = []
        for label_values in range(0, 24, 2):
            xticklabels.append('')
            xticklabels.append(f"{label_values}h")
        self.ax1.set_xticklabels(xticklabels, color='b', alpha=0.6)

    def add_block(self, ra_cen: Quantity, width: Quantity, color='b', alpha=0.1, label=None):
        """
        Add a 'wedge' to the polar plot to show the location of a block.
        Args:
            ra_cen: the Right Ascension of the centre of the block
            width: the full width of the block
        """
        self.ax1.bar(ra_cen.to('rad').value,
                     self.outer_edge.to('au').value,
                     linewidth=0.1,
                     width=width.to('rad').value,
                     bottom=self.inner_edge.to('au').value,
                     zorder=0,
                     color=color,
                     alpha=alpha)
        if label:
            self.ax1.annotate(label, ra_cen.to('rad').value,
                              self.outer_edge.to('au').inner - 3.,
                              size=15, color=color)


    def add_galatic_plane(self):
        # plot exclusion zones due to Galactic plane: RAs indicate where bar starts, rather than its centre angle
        width = math.radians(3 * 15)
        self.ax1.bar(math.radians(4.5 * 15), extent, width=width, color=plot_fanciness.ALMOST_BLACK, linewidth=0, alpha=0.2)
        self.ax1.bar(math.radians(16.5 * 15), extent, width=width, color=plot_fanciness.ALMOST_BLACK, linewidth=0, alpha=0.2)
        self.ax1.annotate('galactic plane', (math.radians(6.9 * 15), extent - 15), size=10, color='k', alpha=0.45)
        self.ax1.annotate('  galactic plane\navoidance zone', (math.radians(16.9 * 15), extent - 12), size=10, color='k', alpha=0.45)

    # again RAs indicate where bar starts, so subtract half the block width from the block centrepoints
    # and approximate blocks as math.radians(7) degrees wide.

    # No optional on these just yet
    def add_planets(self, date='2020/01/01 10:00:00'):
        for planet in [ephem.Saturn(), ephem.Uranus(), ephem.Neptune(), ephem.Pluto()]:
            planet.compute(ephem.date(date))
            fc = 'k'
            if planet.name == 'Pluto':
                alpha = 0.35
                size = 10
                fs = 5
            else:
                alpha = 0.7
                size = 20
                fs = 10
            self.ax1.scatter(planet.ra, planet.sun_distance,
                       marker='o',
                       s=size,
                       facecolor=fc,
                       edgecolor=fc,
                       alpha=alpha)
            if planet.name == 'Saturn':
                self.ax1.annotate(planet.name, (planet.ra + (math.radians(3)), planet.sun_distance - 2), size=fs)
            elif planet.name == 'Uranus':
                self.ax1.annotate(planet.name, (planet.ra - (math.radians(12)), planet.sun_distance + 1), size=fs)
            else:
                self.ax1.annotate(planet.name, (planet.ra - (math.radians(0.5)), planet.sun_distance + 2), size=fs)

            # Neptune's orbit has e = 0.01, so can get away with a circle
            if planet.name == 'Neptune':
                orb = np.arange(0, 2 * np.pi, (2 * np.pi) / 360)
                self.ax1.plot(orb, np.repeat(planet.sun_distance, len(orb)), color='b', linestyle=':', linewidth=0.4, alpha=0.7)

            return
