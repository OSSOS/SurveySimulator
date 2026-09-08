"""
Make the standard 'Rose Plot' for all surveys we have SSim characterization files for.
"""
import ossssim
from ossssim.plotter import RosePlot
import os
from astropy.table import Table
from astropy.time import Time
from matplotlib import pyplot


def face_down_plot() -> None:
    """_
    Plot the detected objects in a face-down plot
    """

    # create a plot object (requires an epoch to put the planets in the right place.)
    plot = RosePlot(Time('2022-09-01T10:00:00').jd)

    # add some KBOs form a 'model file' which is a file with at least the following columns defined
    # a  e inc  peri  node  M  H
    distant_kbos = ossssim.ModelFile('distant_TNOs.txt')
    plot.add_model(distant_kbos, mc='r', ms=10)

    # add pointing wedges from the OSSOS++ surveys
    plot.add_pointings(os.path.join('Surveys/CFEPS/pointings.list'))

    # Add some objects to the plot, this files is in SSim CDS format.
    detection_table = Table.read(os.path.join('Surveys/CFEPS/CFEPS.CDS'), format='cds')
    plot.add_detections(detection_table)

    # put a wedge where the galactic plane is obscuring detections.
    plot.add_galactic_plane()

    # show the plot.. could do pyplot.savefig('xyz.pdf') instead.
    pyplot.savefig('ossos_rose.pdf')


if __name__ == '__main__':
    face_down_plot()
