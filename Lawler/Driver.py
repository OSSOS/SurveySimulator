#! /usr/bin/env python

import sys
from random import random, gauss
import numpy as np
import math


import SurveySubs
import ssimTools as ss
import GiMeObj as go

drad = np.pi/180.0



# Read in driver arguments
with open ('input.file','r') as f:
    tmp = f.read().split()
    distri_file = tmp[1]                  # Name of model file for GiMeObj     
    fsize = tmp[0]                        # Number of rows in file

with open('seeds.txt','r') as f:
    seed = int(f.read())                  # Seed for random number generator

with open('number_to_track.txt','r') as f:
    n_track_max = int(f.readline())           # Number of objects to track

with open('number_to_detect.txt','r') as f:
    n_detect_max = int(f.readline())           # Number of objects to track 

with open ('surveydir.txt','r') as f:
    survey_dir = f.read().split()[0]      # Path to directory containing the characterization files

detect_file = 'detections.dat'            # Output file for detections
track_file = 'tracked.dat'                # Output file for tracked objects

go.setRand(seed)
ss.setRand(seed)

f_detect = ss.detFile(detect_file, seed) # set-up the detection file comments in OSSOS format
f_track = ss.trackFile(track_file) # set-up the tracked file comments in OSSOS format


drawn = open('drawn.dat','w') # the first 5000 objects drawn

drawn.write('# a\t\te\t\tinc\t\tnode\t\tperi\t\tManom\t\tH\t\tresamp \n')

keep_going = True
n_iter, n_hits, n_track = 0, 0, 0

#n_track_max=2

# I believe this is envisioned as a comment for the kind of object returned by GiMeObj
comments = 'res'

while keep_going:

    # Draw an object 
    a,e,inc,node,peri,M,epoch,h,color,gb,ph,period,amp,resamp = go.GiMeObj(distri_file)

    # Write out the first 5000 objects to a file to give a small representative sample
    if n_iter <5000:
        drawn.write(str(a)+'\t'+str(e)+'\t'+str(inc/drad)+'\t'+str(node/drad)+'\t'+str(peri/drad)+'\t'+str(M/drad)+'\t'+str(h)+'\t'+str(resamp)+'\n')

    # Counter: advantage of Python over Fortran: integers can be of any value
    # There is no limit at 2**31-1.
    n_iter += 1

    # Call the survery simulator
    # The output seed2 is never used, but is returned by Fortran so it is stored
    seed2,flag,ra,dec,d_ra,d_dec,r,delta,m_int,m_rand,eff,isur,mt,epochp,ic,surna,h_rand = SurveySubs.detos1(a,e,inc,node,peri,M,epoch,h,color,gb,ph,period,amp,survey_dir,seed)

    # Condition for CFEPS objects with d<20
    if (flag > 0) and ((surna[0]=='L') or (surna[0]=='p')) and (r<20):
        continue


    # If an object is detected, flag > 0
    if flag > 0:
        n_hits += 1
        # Write the detected object out to the detected in the OSSOS format
        ss.detWrite(detect_file, a, e, inc, node, peri, M, resamp, r, mt, m_rand, h_rand, color, ic, flag, delta, m_int, h, eff, ra, dec, d_ra, d_dec, surna, comments)
        # If an object is also tracked
        if (flag > 2) and (np.mod(flag,2) == 0):
            n_track += 1
            # Write the detected object out to the tracked file in the OSSOS format
            ss.trackWrite(track_file, a, e, inc, node, peri, M, resamp, r, mt, m_rand, h_rand, color, ic, comments)


    # break for detections file
    if ((n_detect_max > 0) & (n_hits >= n_detect_max)) | ((n_detect_max < 0) & (n_iter >= -n_track_max)):
        keep_going = False




    # If the break condition of max tracked or max observed is reached
    if ((n_track_max > 0) & (n_track >= n_track_max)) | ((n_track_max < 0) & (n_iter >= -n_track_max)):
        keep_going = False


# Done going through model. Write out summary and quit
ss.detSuffix(detect_file, n_iter, n_hits, n_track)


