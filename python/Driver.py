#! /usr/bin/env python
"""Driver.py is the core of the Survey Simulator's python wrapper. Driver is the module that is executed by the user
which in turn calls all other required modules to run the survey simulator. """

# Module that allows dynamic import capabilities. Using this allows the user to supply arbitrary GiMeObj-type files
# in "input.file" without the need to alter hardcoded values inside Driver.
import importlib

import SurveySubs  # Imports the compiled fortran library.
import numpy as np

import ssimTools as Tools  # Imports functions to format and write output files.

# Read in required parameters from input files.
obj_file, seed, n_detect_max, n_track_max, survey_dir = Tools.input_variables("input.file")

# Prepare Output files.
drawn_file = 'drawn.dat'  # Output file for drawn objects.
detect_file = 'detections.dat'  # Output file for detected objects.
track_file = 'tracked.dat'  # Output file for tracked objects.

Tools.drawn_header(drawn_file)  # Create the tracked file and write the header.
Tools.det_header(detect_file, seed)  # Create the detection file and write the header.
Tools.track_header(track_file)  # Create the tracked file and write the header.

# Import the arbitrarily named GiMeObj module given as the first value in "input.file".
# Name must be given without a file extension.
go = importlib.import_module(obj_file)

# Initialize the seed and distribution arrays in GiMeObj.
go.set_seed(seed)
go.initialize()

n_drawn, n_det, n_track = 0, 0, 0  # Initialize counters.

while (n_det < n_detect_max) and (n_track < n_track_max):

    # Draw an object 
    a, e, inc, node, peri, M, epoch, h, color, gb, ph, period, amp, resamp, comments = go.gimeobj()

    # Write out the first 5000 objects to a file to give a small representative sample
    if n_drawn < 5000:
        Tools.drawn_write(
            drawn_file, a, e, inc, node, peri, M, h, resamp, comments
        )

    n_drawn += 1

    # Call the survey simulator
    # The output seed2 is never used, but is returned by Fortran so it is stored
    seed2, flag, ra, dec, d_ra, d_dec, r, delta, m_int, m_rand, eff, isur, mt, epochp, ic, surna, h_rand = \
        SurveySubs.detos1(a, e, inc, node, peri, M, epoch, h, color, gb, ph, period, amp, survey_dir, seed)

    # Condition for CFEPS objects with d<20. The survey simulator doesn't return objects that would not have been
    # seen by OSSOS survey blocks. CFEPS has a different limit that isn't checked by the survey simulator and must be
    # filtered for here. Since such an object would not have been detected by the CFEPS survey the detected and
    # tracked section of Driver are skipped when this condition is true.
    if (flag > 0) and ((surna[0] == 'L') or (surna[0] == 'p')) and (r < 20):
        continue

    # If an object is detected, flag > 0
    if flag > 0:
        n_det += 1
        # Write the detected object out to the detected file.
        Tools.det_write(
            detect_file, a, e, inc, node, peri, M, resamp, r, mt, m_rand, h_rand, color,
            ic, flag, delta, m_int, h, eff, ra, dec, d_ra, d_dec, surna, comments
        )
        # If an object is also tracked
        if (flag > 2) and (np.mod(flag, 2) == 0):
            n_track += 1
            # Write the detected object out to the tracked file.
            Tools.track_write(
                track_file, a, e, inc, node, peri, M, resamp, r, mt, m_rand, h_rand, color, ic, comments
            )

# Done going through model. Write out summary and quit
Tools.det_suffix(detect_file, n_drawn, n_det, n_track)
