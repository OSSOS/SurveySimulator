"""
Example of how to make some basic plots using the ossssim.plotting module
"""
import ossssim
from ossssim import plotter

# read in the L7 model as the background model for plot
model = ossssim.ModelFile('../../Models/L7model-3.0-9.0')
# read in the detection file to show locations of model detections.
detect = ossssim.ModelFile('L7detect.dat')

# Create a FaceDownPlot instance
plot_driver = plotter.FaceDownPlot(longitude_neptune=model.longitude_neptune)
# plot some scale rings
plot_driver.plot_rings()
# plot the background Model, plots 10**4 model randomly selected objects.
plot_driver.add_model(model, sample_size=10**4)
# plot the detected objects, plots all objects but in red with larger points
plot_driver.add_model(detect, mc='r', ms=5.)
# show the plot (could do plt.savefig here instead)
plot_driver.show()
