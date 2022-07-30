"""
Make the standard 'Rose Plot' for all surveys we have SSim characterization files for.
"""
import ossssim
import os
from astropy.table import Table
from astropy.time import Time

SURVEYS_DIR = '/arc/projects/classy/Surveys'


def face_down_plot() -> None:
    """_
    Plot the detected objects in a face-down plot
    """

    # create a plot object (requires an epoch to put the planets in the right place.)
    plot = ossssim.RosePlot(Time('2022-09-01T10:00:00').jd)

    # add some KBOs form a 'model file' which is a file with at least the following columns defined
    # a  e inc  peri  node  M  H
    distant_kbos = ossssim.ModelFile('distant_TNOs.txt')
    plot.add_model(distant_kbos, mc='r', ms=10)

    # add pointing wedges from the OSSOS++ surveys
    plot.add_pointings(os.path.join(SURVEYS_DIR, 'ObsSummary/All_Surveys/pointings.list'))

    # Add the classy pointing plan.
    classy_pointings_list = os.path.join(SURVEYS_DIR, 'ObsSummary/CLASSY/pointings.list')
    plot.add_pointings(classy_pointings_list, color='g', alpha=0.7, label=True)

    # Add the Deep pointing plan.
    deep_pointings_list = os.path.join(SURVEYS_DIR, 'ObsSummary/DEEP/pointings.list')
    plot.add_pointings(deep_pointings_list, color='r', alpha=0.3)

    # Add some objects to the plot, this files is in SSim CDS format.
    detection_table = Table.read(os.path.join(SURVEYS_DIR,
                                              'ObsSummary/All_Surveys/All_Surveys_v11.CDS'),
                                 format='cds')
    plot.add_detections(detection_table)

    # put a wedge where the galactic plane is obscuring detections.
    plot.add_galactic_plane()

    # show the plot.. could do pyplot.savefig('xyz.pdf') instead.
    plot.show()


if __name__ == '__main__':
    face_down_plot()
