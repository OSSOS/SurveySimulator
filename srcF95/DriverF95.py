#! /usr/bin/env python3

import sys, time
sys.path.append('.')

def conv(st):
    return "".join(chr(x) for x in st)

import numpy as num
drad = num.pi/180.0

import SurveySubsF95
import GiMeObjF95

SS = SurveySubsF95.surveysub
GO = GiMeObjF95.gimeobjut
color = num.zeros(10)

# Read in driver arguments
# Seed for random number generator
stri = input()
seed = int(stri.split()[0])
# Maximum number of detections (>0) or -maximum number of trials (<0)
stri=input()
n_track_max = int(stri.split()[0])
# Directory containing the characterization files
survey_dir = input()
# File with model parameters
distri_file = input()
# Output file for detections
detect_file = input()
# Output file for tracked objects
track_file = input()

# Open detection file and write header
f_detect = open(detect_file, 'w')
f_detect.write ('# Seed: %10d\n#\n' % (seed))
dati = time.strftime("%Y-%m-%dT%H:%M:%S.000  %z")
f_detect.write ('# Creation time: {:30s}\n#\n'.format(dati))
f_detect.write ('# flag: >0: detected; >2: characterized; 0 mod(2): tracked\n')
f_detect.write ('# Survey: name of the block\n#\n')
f_detect.write ('#   a      e        i        q        r        M       node     peri  m_rand H_rand color flag delta    m_int    H_int eff   RA(H)     DEC    Surv. Comments\n')

# Open tracked detection file (no header)
f_track = open(track_file, 'w')
f_track.write ('# Creation time: {:30s}\n#\n'.format(dati))
f_track.write ('#   a      e        i        q        r        M       node     peri  m_rand H_rand color Comments\n')

keep_going = True
n_iter, n_hits, n_track = 0, 0, 0

while keep_going:
    nchar = 0
    o_m, epoch, h, gb, ph, period, amp, commen, nchar, ierr = GO.gimeobj(distri_file, seed, color)
    if ierr == -20:
        print('GiMeObj returned -20, stopping.')
        break
    if ierr == 100:
        keep_going = False
    if ierr != -10:
# Counter: advantage of Python over Fortran: integers can be of any value
# There is no limit at 2**31-1.
        n_iter += 1

        flag, ra, dec, d_ra, d_dec, r, delta, m_int, m_rand, eff, isur, mt, jdayp, ic, surna, h_rand = SS.detos1(o_m, epoch, h, color, gb, ph, period, amp, survey_dir, seed)

        if flag > 0:
# m_int and h are in "x" band (filter of object creation)
# m_rand and h_rand are in discovery filter
            n_hits += 1
            f_detect.write('{:8.3f} {:6.3f} {:8.3f} {:8.3f} {:8.3f} {:8.3f} {:8.3f} {:8.3f} {:6.2f} {:6.2f} {:5.2f} {:2d} {:8.3f} {:8.3f} {:6.2f} {:4.2f} {:8.5f} {:8.4f} {:6s} {:s}\n'.format(o_m.a, o_m.e, o_m.inc/drad, o_m.a*(1.-o_m.e), r, mt/drad, o_m.node/drad, o_m.peri/drad, m_rand, h_rand, color[ic-1], flag, delta, m_int, h, eff, ra/drad/15., dec/drad, conv(surna), conv(commen[0:nchar])))
            if (flag > 2) and (num.mod(flag,2) == 0):
                n_track += 1
                f_track.write("{:8.3f} {:6.3f} {:8.3f} {:8.3f} {:8.3f} {:8.3f} {:8.3f} {:8.3f} {:6.2f} {:6.2f} {:5.2f} {:s}\n".format(o_m.a, o_m.e, o_m.inc/drad, o_m.a*(1.-o_m.e), r, mt/drad, o_m.node/drad, o_m.peri/drad, m_rand, h_rand, color[ic-1], conv(commen[0:nchar])))

    if ((n_track_max > 0) & (n_track >= n_track_max)) | ((n_track_max < 0) & (n_iter >= -n_track_max)):
        keep_going = False

# Done going through model. Write out summary and quit
f_detect.write ('# Total number of objects:   {:11d}\n'.format(n_iter))
f_detect.write ('# Number of detections:      {:7d}\n'.format(n_hits))
f_detect.write ('# Number of tracked objects: {:7d}\n'.format(n_track))
f_detect.close()
f_track.close()
