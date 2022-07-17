"""
Example of how to make some basic plots using the osssim.plotting module
"""
import osssim
from osssim import plotter
import numpy

# first we run a simulation, so we have something interesting to plot.
model = osssim.SSimModelFile('../../../Models/L7model-3.0-9.0')
ssim = osssim.OSSSSim('../../../Surveys/CFEPS')

seed = 113715191
track_condition = 0

# create a detection output file and fill header with circumstance values that match the input model.
"""
detect_file = osssim.DetectFile('L7detect.dat')
detect_file.epoch = model.epoch
detect_file.lambda_neptune = model.longitude_neptune
detect_file.colors = model.colors
detect_file.write_header(seed)

track_file = osssim.DetectFile('L7track.dat')
track_file.epoch = model.epoch
track_file.lambda_neptune = model.longitude_neptune
track_file.colors = model.colors
track_file.write_header(seed)

# setup some variables to count the number of object tracked.
n_iter = n_hits = n_track = 0

# set track condition to stop when we have 1 tracked object
# set value to -1 to run until 1 iterations has occurred
# set value to 0 to loop to the end of the file.
for row in model:
    n_iter += 1
    result = ssim.simulate(row, seed=seed, epoch=model.epoch)
    if result['flag'] > 0:
        n_hits += 1
        detect_file.write_row(result)

    if (result['flag'] > 2) and (numpy.mod(result['flag'], 2) == 0):
        n_track += 1
        track_file.write_row(result)
    # stop the loop when maximum/desired detections have occurred if ntrack==0 then go to end of model file.
    if (0 < track_condition <= n_track) | (0 < -track_condition <= n_iter):
        break
detect_file.write_footer(n_iter=n_iter, n_hits=n_hits, n_track=n_track)
track_file.write_footer(n_iter=n_iter, n_hits=n_hits, n_track=n_track)
"""
# read back in the detection file as a Model File as that's what the plotter wants.
detect = osssim.SSimModelFile('L7detect.dat')
plot_driver = plotter.FacedownPlot(longitude_neptune=model.longitude_neptune)
plot_driver.plot_rings()
plot_driver.add_model(model, mc='k', sample_size=10**4)
plot_driver.add_model(detect, mc='r', ms=5.)
plot_driver.show()