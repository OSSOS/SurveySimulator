"""
Some methods to aid making orbit plots from simulator ouptuts
"""
import osssim
from matplotlib import pyplot as plt
import numpy
from numpy.random import default_rng
from . import definitions
import logging
np = numpy

class TimeSeriesPlot(object):

    def __init__(self, model):
        """Given an input model setup methods to make plots.
        Args:
            model (osssim.SSimModelFile or osssim.SSimDetect): a model object from the osssim package."""
        self.model = model
        self.fig = plt.figure(figsize=(8,15))

    def plot(self, variables=None, samples=10**6):
        """
        Plot model elements as the value vs position in the array (index).

        Args:
            variables ([] or None): a list of columns to plot. e.g. (['a', 'e', 'inc']) default:None plots all columns in model.
            samples (int): upto samples points will be plotted.

        To see what columns are available look at self.model.colnames

        Intended for diagnostic purposes.
        """

        if variables is None:
            variables = self.model.targets.colnames
        nx = len(variables)//2
        ny = 2
        iter = 0
        for column in variables:
            iter += 1
            ax = self.fig.add_subplot(nx,ny,iter)
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


class FacedownPlot(object):

    def __init__(self, longitude_neptune):
        """Given an input model setup methods to make plots.
        Args:
            model (osssim.SSimModelFile or osssim.SSimDetect): a model object from the osssim package."""
        self.fig = plt.figure(figsize=(8, 8))
        self.ax = self.fig.add_subplot(111)
        self.lambdaN = longitude_neptune

    def plot_rings(self):

        # Neptune's position
        xN = 30.1 * np.cos(self.lambdaN)
        yN = 30.1 * np.sin(self.lambdaN)

        self.ax.plot(xN, yN, 'ob', mec='k')
        self.ax.plot(0, 0, 'oy', mec='k')
        for guide_circle in [10., 30., 50., 100.]:
            circle = plt.Circle((0, 0), guide_circle, fill=False, ec='b', ls=':')
            self.ax.add_artist(circle)
            self.ax.text(0, guide_circle+1, f"{guide_circle} au",
                     horizontalalignment='center', color='b')

    def add_model(self, model, mc='k', ms=1., sample_size=None):
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

        self.ax.plot(x, y, f'.{mc}', markersize=ms)

    def show(self):
        plt.ylim(-100,100)
        plt.xlim(-100,100)
        plt.axis('off')
        plt.show()
