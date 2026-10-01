"""
Example: basic RosePlot of a model file.

Lookup-table models (e.g. L7) provide Keplerian elements only; RosePlot
converts those to heliocentric x,y,z for plotting.
"""
from pathlib import Path
from matplotlib import pyplot as plt
import ossssim
from ossssim import plotter
from astropy.time import Time

REPO = Path(__file__).resolve().parents[1]


def run():
    model_path = REPO / 'F95' / 'tests' / 'Models' / 'L7model-3.0-9.0'
    if not model_path.is_file():
        # Fall back to sibling SurveySimulator-Data if present
        model_path = REPO.parent / 'SurveySimulator-Data' / 'Models' / 'L7model-3.0-9.0'
    model = ossssim.ModelFile(str(model_path))

    plot_driver = plotter.RosePlot(epoch=Time(model.epoch))
    plot_driver.add_scale_rings()
    plot_driver.add_model(model, sample_size=10**4)
    # plot_driver.show()
    plt.savefig('roseplot.png')


if __name__ == '__main__':
    run()
