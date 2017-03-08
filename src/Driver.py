#! /usr/bin/env python

import sys
sys.path.append('.')

import numpy as num
drad = num.pi/180.0

import SurveySubs
import GiMeObj

# Read in driver arguments
# Seed for random number generator
str = raw_input()
seed = int(str.split()[0])
# Maximum number of detections (>0) or -maximum number of trials (<0)
str=raw_input()
n_track_max = int(str.split()[0])
# Directory containing the characterization files
survey_dir = raw_input()
# File with model parameters
distri_file = raw_input()
# Output file for detections
detect_file = raw_input()
# Output file for tracked objects
track_file = raw_input()

# Open detection file and write header
f_detect = open(detect_file, 'w')
f_detect.write ('# Seed: %10d\n#\n' % (seed))
f_detect.write ('# flag: >0: detected; >2: characterized; 0 mod(2): tracked\n')
f_detect.write ('# Survey: name of the block\n')
f_detect.write ('# delta_ra: distance from center of pointing [arcsec]\n')
f_detect.write ('# delt_dec: distance from center of pointing [arcsec]\n#\n')
f_detect.write ('#   a      e        i        q        r        M      m_rand H_rand color flag delta    m_int   H_int eff   RA(H)     DEC    delta_ra delt_dec Surv.  Comments\n')

# Open tracked detection file (no header)
f_track = open(track_file, 'w')
f_track.write ('#   a      e        i        q        r        M      m_rand H_rand color Comment\n')

keep_going = True
n_iter, n_hits, n_track = 0, 0, 0

while keep_going:
    nchar = 0
    seed,a,e,inc,node,peri,M,epoch,h,color,gb,ph,period,amp,comments,nchar,ierr = GiMeObj.gimeobj(distri_file,seed)
    if ierr == -10:
        pass
    if ierr == -20:
        break
    if ierr == 100:
        keep_going = False

# Counter: advantage of Python over Fortran: integers can be of any value
# There is no limit at 2**31-1.
    n_iter += 1

    seed,flag,ra,dec,d_ra,d_dec,r,delta,m_int,m_rand,eff,isur,mt,epochp,ic,surna,h_rand = SurveySubs.detos1(a,e,inc,node,peri,M,epoch,h,color,gb,ph,period,amp,survey_dir,seed)

    if flag > 0:
# m_int and h are in "x" band (filter of object creation)
# m_rand and h_rand are in discovery filter
        n_hits += 1
        f_detect.write('%8.3f %6.3f %8.3f %8.3f %8.3f %8.3f %8.3f %6.2f %5.2f %2d %8.3f %8.3f %6.2f %4.2f %8.5f %8.4f %8.5f %8.5f %6s %s\n' % (a, e, inc/drad, a*(1.-e), r, mt/drad, m_rand, h_rand, color[ic-1], flag, delta, m_int, h, eff, ra/drad/15., dec/drad, d_ra/drad*3600./24., d_dec/drad*3600./24., surna, comments[0:nchar]))
        if (flag > 2) and (num.mod(flag,2) == 0):
            n_track += 1
            f_track.write("%8.3f %6.3f %8.3f %8.3f %8.3f %8.3f %8.3f %6.2f %5.2f %s\n" % (a, e, inc/drad, a*(1.-e), r, mt/drad, m_rand, h_rand, color[ic-1], comments[0:nchar]))

    if ((n_track_max > 0) & (n_track >= n_track_max)) | ((n_track_max < 0) & (n_iter >= -n_track_max)):
        keep_going = False

# Done going through model. Write out summary and quit
f_detect.write ('# Total number of objects:   %11d\n' % (n_iter))
f_detect.write ('# Number of detections:      %7d\n' % (n_hits))
f_detect.write ('# Number of tracked objects: %7d\n' % (n_track))
f_detect.close()
f_track.close()
