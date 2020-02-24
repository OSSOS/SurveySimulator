#! /usr/bin/env python

import importlib

import SurveySubs
import numpy as np

import ssimTools as Tools

drad = np.pi / 180.0

# Read in driver arguments
with open('input.file', 'r') as f:
    tmp = f.read().split()
    distri_file = tmp[1]  # Name of model file for GiMeObj
    obj_file = tmp[0]  # Name of GiMeObj module file (without .py file extension)

with open('seeds.txt', 'r') as f:
    seed = int(f.read())  # Seed for random number generator

with open('number_to_track.txt', 'r') as f:
    n_track_max = int(f.readline())  # Number of objects to track

with open('number_to_detect.txt', 'r') as f:
    n_detect_max = int(f.readline())  # Number of objects to track

with open('surveydir.txt', 'r') as f:
    survey_dir = f.read().split()[0]  # Path to directory containing the characterization files

detect_file = 'detections.dat'  # Output file for detections
track_file = 'tracked.dat'  # Output file for tracked objects

# Import the arbitrarily named GiMeObj module given in input.file. Name must be given without a file extension.
go = importlib.import_module(obj_file)

go.setrand(seed)
Tools.setrand(seed)

f_detect = Tools.detfile(detect_file, seed)  # set-up the detection file comments in OSSOS format
f_track = Tools.trackfile(track_file)  # set-up the tracked file comments in OSSOS format

drawn = open('drawn.dat', 'w')  # the first 5000 objects drawn

drawn.write(f'#{"a":>7} {"e":>6} {"inc":>8} {"node":>8} {"peri":>8} {"Manom":>8} {"H":>6} {"resamp":>8}\n')

keep_going = True
n_iter, n_hits, n_track = 0, 0, 0

# n_track_max=2

# I believe this is envisioned as a comment for the kind of object returned by GiMeObj
comments = 'res'

while keep_going:

    # Draw an object 
    a, e, inc, node, peri, M, epoch, h, color, gb, ph, period, amp, resamp = go.gimeobj(distri_file)

    # Write out the first 5000 objects to a file to give a small representative sample
    if n_iter < 5000:
        drawn.write(
            f"{a:8.3f} {e:6.3f} {inc/drad:8.3f} {node/drad:8.1f} {peri/drad:8.1f} {M/drad:8.1f} {h:6.2f} {resamp:8.1f}"
            f"\n"
        )

    # Counter: advantage of Python over Fortran: integers can be of any value
    # There is no limit at 2**31-1.
    n_iter += 1

    # Call the survery simulator
    # The output seed2 is never used, but is returned by Fortran so it is stored
    seed2, flag, ra, dec, d_ra, d_dec, r, delta, m_int, m_rand, eff, isur, mt, epochp, ic, surna, h_rand = \
        SurveySubs.detos1(a, e, inc, node, peri, M, epoch, h, color, gb, ph, period, amp, survey_dir, seed)

    # Condition for CFEPS objects with d<20
    if (flag > 0) and ((surna[0] == 'L') or (surna[0] == 'p')) and (r < 20):
        continue

    # If an object is detected, flag > 0
    if flag > 0:
        n_hits += 1
        # Write the detected object out to the detected in the OSSOS format
        Tools.detwrite(
            detect_file, a, e, inc, node, peri, M, resamp, r, mt, m_rand, h_rand, color,
            ic, flag, delta, m_int, h, eff, ra, dec, d_ra, d_dec, surna, comments
        )
        # If an object is also tracked
        if (flag > 2) and (np.mod(flag, 2) == 0):
            n_track += 1
            # Write the detected object out to the tracked file in the OSSOS format
            Tools.trackwrite(track_file, a, e, inc, node, peri, M, resamp, r, mt, m_rand, h_rand, color, ic, comments)

    # break for detections file
    if ((n_detect_max > 0) & (n_hits >= n_detect_max)) | ((n_detect_max < 0) & (n_iter >= -n_track_max)):
        keep_going = False

    # If the break condition of max tracked or max observed is reached
    if ((n_track_max > 0) & (n_track >= n_track_max)) | ((n_track_max < 0) & (n_iter >= -n_track_max)):
        keep_going = False

# Done going through model. Write out summary and quit
Tools.detsuffix(detect_file, n_iter, n_hits, n_track)
