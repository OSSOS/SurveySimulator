"""
Some methods to aid making orbit plots from simulator outputs
"""
import logging
import numpy
from matplotlib import pyplot as plt
from numpy.random import default_rng

from . import definitions
from .models import ModelFile
from .models import Parametric


np = numpy


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


class FaceDownPlot:
    """
    Plot objects of model as face-down view of the outer solar system.
    """
    def __init__(self, longitude_neptune: float) -> None:
        """Given an input model setup methods to make plots.
        Args:
            longitude_neptune: Neptune's longitude at the epoch of the plot.
        """
        self.fig = plt.figure(figsize=(8, 8))
        self.ax = self.fig.add_subplot(111)
        self.lambdaN = longitude_neptune

    def plot_rings(self, radii: list = (10, 30, 50, 100)) -> None:
        """
        Add scale rings to the face-down plot.

        radii: the radii of rings to plot to provide scale to the plot.
        """

        # Neptune's position
        xN = 30.1 * np.cos(self.lambdaN)
        yN = 30.1 * np.sin(self.lambdaN)

        self.ax.plot(xN, yN, 'ob', mec='k')
        self.ax.plot(0, 0, 'oy', mec='k')
        for guide_circle in radii:
            circle = plt.Circle((0, 0), guide_circle, fill=False, ec='b', ls=':')
            self.ax.add_artist(circle)
            self.ax.text(0, guide_circle + 1, f"{guide_circle} au",
                         horizontalalignment='center', color='b')

    def add_model(self, model: ModelFile, mc: str = 'k', ms: float = 1., sample_size: int = None, alpha: float = 1.0) -> None:
        """
        Make a face-down plot of the solar system for this model.
        """
        state = model.cartesian
        if sample_size is None:
            sample_size = len(model.targets)
        rng = default_rng()
        choice = rng.integers(0, len(model.targets), sample_size)
        x = state['x'][choice]
        y = state['y'][choice]
        # z = state['z'][choice]

        self.ax.plot(x, y, f'.{mc}', markersize=ms, alpha=alpha)

    def show(self) -> None:
        """
        Display the plot.
        """
        self.ax.set_ylim(-100, 100)
        self.ax.set_xlim(-100, 100)
        plt.axis('off')
        plt.show()
