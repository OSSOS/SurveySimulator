"""
Helpers for orbit / survey footprint plots from simulator outputs.
"""
from __future__ import annotations

import logging
from pathlib import Path
from typing import List, Optional, Union

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
from .pos_cart import pos_cart
from .survey import SurveyCharacterization

np = numpy

# setup the plotting Look and Feel.
rcParams['font.size'] = 12  # good for posters/slides
rcParams['patch.facecolor'] = (0.4, 0.7607843137254902, 0.6470588235294118)
rcParams['figure.figsize'] = (10, 10)
rcParams['figure.dpi'] = 150
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
# Prefer common fonts; fall back silently if Sofia Pro is missing
rcParams['font.sans-serif'] = ['DejaVu Sans', 'Tahoma', 'Helvetica', 'Arial', 'Sofia Pro']
ALMOST_BLACK = '#262626'


def _as_time(epoch) -> Time:
    """Normalize constructor epoch to an astropy Time."""
    if isinstance(epoch, Time):
        return epoch
    if isinstance(epoch, (int, float, numpy.floating)):
        return Time(epoch, format='jd')
    if isinstance(epoch, Quantity):
        raise TypeError(
            'RosePlot epoch must be an astropy Time (or JD float), '
            'not an angle Quantity such as Neptune longitude'
        )
    return Time(epoch)


# Angle (degrees) along which radial distance labels are drawn, and a
# semi-transparent white box so labels stay legible over dense scatter.
_RADIAL_LABEL_ANGLE = 100
_RADIAL_LABEL_BBOX = dict(boxstyle='round,pad=0.15', facecolor='white',
                          alpha=0.6, edgecolor='none')


def _wedge_width_from_area(area_deg2: float) -> Quantity:
    """Approximate RA wedge width from spherical area (square-equivalent)."""
    side = float(numpy.sqrt(max(area_deg2, 0.0)))
    return side * units.deg


def _ensure_cartesian_table(table):
    """
    Return a table that has heliocentric ecliptic ``x,y,z``.

    Detect/parametric outputs already store Rebound state vectors. Lookup-table
    models (e.g. L7) only have Keplerian elements; convert with the Python
    ``pos_cart`` port of Fortran ``elemutils.pos_cart`` (not currently
    f90wrap-exported — only datadec/ioutils/surveysub are wrapped).
    """
    if all(name in table.colnames for name in ('x', 'y', 'z')):
        return table
    required = ('a', 'e', 'inc', 'node', 'peri', 'M')
    missing = [name for name in required if name not in table.colnames]
    if missing:
        raise KeyError(
            "RosePlot.add_model needs cartesian columns x,y,z or Keplerian "
            f"elements {required}; missing {missing}"
        )

    def _col(name, unit=None):
        col = table[name]
        if unit is not None and hasattr(col, 'to'):
            return numpy.asarray(col.to(unit).value, dtype=float)
        if hasattr(col, 'value'):
            return numpy.asarray(col.value, dtype=float)
        return numpy.asarray(col, dtype=float)

    xyz = pos_cart(
        _col('a', 'au'),
        _col('e'),
        _col('inc', 'rad'),
        _col('node', 'rad'),
        _col('peri', 'rad'),
        _col('M', 'rad'),
    )
    out = table.copy()
    out['x'] = xyz[0] * units.au
    out['y'] = xyz[1] * units.au
    out['z'] = xyz[2] * units.au
    return out


class TimeSeriesPlot:
    """
    Plot model elements vs index (diagnostic distribution checks).
    """

    def __init__(self, model: (ModelFile or Parametric)) -> None:
        self.model = model
        self.fig = plt.figure(figsize=(8, 15))

    def plot(self, variables: list = None) -> None:
        if variables is None:
            variables = self.model.targets.column_names
        n = max(len(variables), 1)
        nx = max((n + 1) // 2, 1)
        ny = 2
        for i, column in enumerate(variables):
            ax = self.fig.add_subplot(nx, ny, i + 1)
            if column not in self.model.targets.column_names:
                logging.warning(f"Could not plot {column} as does not appear input model.")
                continue

            values = self.model.targets[column]
            if isinstance(values[0], (list, numpy.ndarray)):
                for idx in range(len(values[0])):
                    ax.plot(values[:, idx].to(definitions.column_unit[column]).value,
                            color='k', marker='o', linestyle='none', linewidth=2, markersize=1)
            else:
                ax.plot(values.to(definitions.column_unit[column]).value,
                        color='k', marker='o', linestyle='none', linewidth=2, markersize=1)
            ax.set_ylabel(f"{column} ({definitions.column_unit[column]})")
        plt.show()


class RosePlot:
    """
    Face-down (polar RA × heliocentric distance) view of models, detections,
    and survey footprints.

    Survey footprints are loaded via ``SurveyCharacterization`` (Fortran
    ``survey_load`` / pointing geometry), not a custom ``pointings.list`` parser.
    """

    def __init__(self, epoch, outer_edge=85 * units.au, inner_edge=10 * units.au) -> None:
        """
        Args:
            epoch: discovery / plot epoch as ``astropy.time.Time`` or Julian Date float.
        """
        self.epoch = _as_time(epoch)
        self._longitude_neptune = None
        self.frame = 'heliocentrictrueecliptic'
        self.outer_edge = outer_edge.to('au').value
        self.inner_edge = inner_edge.to('au').value
        self.fig = plt.figure(figsize=(8, 8))
        rect = [0.0725, 0.0725, 0.85, 0.85]
        self.ax1 = self.fig.add_axes(rect, polar=True, frameon=False)
        self.ax1.set_aspect('equal')

        self.ax1.set_rlim(0, self.outer_edge)
        rings = list(range(0, 100, 15))
        ring_labels = [""] * 2 + [f"{x:3d} au" for x in rings[2:]]
        # pad / trim so label count matches rings
        if len(ring_labels) < len(rings):
            ring_labels = ring_labels + [f"{x:3d} au" for x in rings[len(ring_labels):]]
        ring_labels = ring_labels[:len(rings)]
        self.ax1.set_rgrids(rings, ring_labels, angle=_RADIAL_LABEL_ANGLE, alpha=0.45)
        self.ax1.yaxis.set_major_locator(MultipleLocator(25))
        for tick_label in self.ax1.yaxis.get_ticklabels():
            tick_label.set_bbox(_RADIAL_LABEL_BBOX)
        self.ax1.xaxis.set_major_locator(MultipleLocator(numpy.deg2rad(15)))
        self.ax1.grid(axis='x', color='k', linestyle='--', alpha=0.2)
        x_tick_labels = []
        lon = numpy.arange(0, 360, 30) * units.deg
        lat = numpy.zeros(len(lon)) * units.deg
        dist = numpy.ones(len(lon)) * 45 * units.au
        coord = SkyCoord(lon, lat, distance=dist, obstime='2000-01-01',
                         frame='heliocentrictrueecliptic').transform_to('icrs')

        for label_values in coord.ra.hour:
            lv = int(numpy.round(label_values))
            x_tick_labels.append('')
            x_tick_labels.append(f"{lv}h")
        self.ax1.set_xticklabels(x_tick_labels, color='b', alpha=0.6)

    @property
    def longitude_neptune(self):
        """Ecliptic longitude of Neptune at the plot epoch (via Horizons)."""
        if self._longitude_neptune is None:
            planet = jplhorizons.Horizons(899, epochs=self.epoch.jd, location='568')
            eph = planet.ephemerides()
            self._longitude_neptune = eph['EclLon'][0]
        return self._longitude_neptune

    def add_pointings(
        self,
        survey_directory: Union[str, Path],
        color: str = 'b',
        high_latitude_color: str = 'y',
        alpha: float = 0.1,
        label: bool = False,
        latitude_cut: float = 10.0,
    ) -> int:
        """
        Add approximate RA wedges for each pointing in a survey characterization.

        Args:
            survey_directory: path to a characterization directory (``pointings.list``
                + ``.eff`` files), loaded via Fortran ``SurveyCharacterization``.
            color: wedge face color near the ecliptic
            high_latitude_color: color when |ecliptic lat| exceeds ``latitude_cut``
            alpha: wedge transparency
            label: annotate with short pointing id
            latitude_cut: degrees; above this |b| use ``high_latitude_color``

        Returns:
            Number of wedges drawn.
        """
        survey = SurveyCharacterization.from_directory(str(survey_directory))
        names: List[str] = []
        n_drawn = 0
        for pointing in survey:
            pos = SkyCoord(
                ra=pointing.ra * units.rad,
                dec=pointing.dec * units.rad,
                distance=44 * units.au,
                obstime='2000-01-01',
            ).transform_to(self.frame)
            wedge_color = color
            if not (-latitude_cut < pos.lat.degree < latitude_cut):
                wedge_color = high_latitude_color
            width = _wedge_width_from_area(pointing.area_deg2)
            name = None
            if label:
                name = pointing.id[:8]
                if name in names:
                    name = None
                else:
                    names.append(name)
            self.add_block(pos.lon, width, color=wedge_color, alpha=alpha, label=name)
            n_drawn += 1
        return n_drawn

    def add_block(self, ra_cen: Quantity, width: Quantity, color='b',
                  alpha=0.1, label=None) -> None:
        """Add a polar wedge for one survey block."""
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
        """Shade ecliptic longitudes where |b| < ``minimum_latitude``."""
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

    def add_planets(self):
        """Plot major planets at the plot epoch via Horizons."""
        ids = {'Jupiter': 599, 'Saturn': 699, 'Uranus': 799,
               'Neptune': 899, 'Pluto': 999}
        fc = ALMOST_BLACK
        for planet_name in ['Jupiter', 'Saturn', 'Uranus', 'Neptune', 'Pluto']:
            planet = jplhorizons.Horizons(ids[planet_name],
                                          location='568',
                                          epochs=self.epoch.jd)
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
        """Plot CDS-style detections (RAdeg, DEdeg, dist)."""
        coords = SkyCoord(detection_table['RAdeg'],
                          detection_table['DEdeg'],
                          obstime='2000-01-01',
                          distance=detection_table['dist'],
                          frame='icrs').transform_to(self.frame)
        self.ax1.scatter(coords.lon.to('rad').value,
                         coords.distance.to('au').value,
                         s=5,
                         c='c')

    def add_scale_rings(self, radii: Optional[List] = None) -> None:
        """Add guide circles at the given heliocentric distances (au)."""
        if radii is None:
            radii = [10, 30, 50, 100]
        theta = numpy.arange(0, 2*numpy.pi, 2*numpy.pi/1000)
        label_angle = numpy.deg2rad(_RADIAL_LABEL_ANGLE)
        for guide_circle in radii:
            r = numpy.ones(len(theta))*guide_circle
            self.ax1.plot(theta, r, ls=':')
            # Label on the ring itself (like the r-grid labels), boxed so it
            # stays readable on top of dense model scatter.
            self.ax1.annotate(f"{guide_circle} au", (label_angle, guide_circle),
                              color='b', horizontalalignment='center',
                              verticalalignment='center',
                              bbox=_RADIAL_LABEL_BBOX)

    def add_model(self, model: ModelFile, mc: str = 'k', ms: float = 1.,
                  sample_size: int = None, alpha: float = 1.0) -> None:
        """
        Scatter model objects on the face-down plot.

        Uses ``x,y,z`` when present (detect / parametric outputs). For
        element-only model files, fills those columns via Keplerian ``pos_cart``.
        """
        table = model.table
        n = len(table)
        if sample_size is not None and sample_size < n:
            rng = default_rng()
            table = table[rng.integers(0, n, sample_size)]
        table = _ensure_cartesian_table(table)

        coord = SkyCoord(table['x'], table['y'], table['z'],
                         representation_type='cartesian',
                         frame='heliocentrictrueecliptic', obstime='2000-01-01').transform_to(self.frame)

        self.ax1.plot(coord.lon.to('rad').value,
                      coord.distance.to('au').value,
                      f'.{mc}',
                      markersize=ms,
                      alpha=alpha)

    @staticmethod
    def savefig(filename, **kwargs) -> None:
        plt.savefig(filename, **kwargs)

    @staticmethod
    def show() -> None:
        plt.show()
