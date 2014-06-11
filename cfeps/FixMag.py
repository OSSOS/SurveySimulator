#! /usr/bin/env python

import sys, os
#! /usr/bin/env python

import sys
sys.path.append('/home/petit/Observations/data/OSSOS')
import SurveySubs as SS

from math import sqrt, tanh, exp, pi
from numpy import argsort

rates = {'o3e' : [0.5, 8.0, 11.0, 15.0],
         'o3o' : [0.5, 7.0, 10.0, 15.0]}
thresholds = {'o3e' : [24.04, 23.96, 23.84],
              'o3o' : [24.40, 24.36, 24.19]}
effs = {'o3e' : [[0.893732548, 2.81616841E-02, 24.1415253, 0.154886752],
                 [0.899760783, 2.31448822E-02, 24.0047779, 0.144415855],
                 [0.878224790, 2.53052488E-02, 23.8959866, 0.152918205]],
        'o3o' : [[0.841142952, 2.08144821E-02, 24.5515594, 0.109849975],
                 [0.881142378, 1.95594411E-02, 24.4109688, 0.115364030],
                 [0.862005889, 1.86797809E-02, 24.2493134, 0.142319098]]}
version = 'OSSOSv3'
ast_dir = '/home/dbase/TNOdb/dbase/data/ast/'
orb_dir = '/home/dbase/TNOdb/dbase/data/orbb/'
n1, n2, n3, nobj = {'o3e' : 0, 'o3o' : 0}, {'o3e' : 0, 'o3o' : 0}, {'o3e' : 0, 'o3o' : 0}, {'o3e' : 0, 'o3o' : 0}
code = 500
gb = -0.12
tol = 1./10. # 6 arcmin
radd = 180./pi
drad = pi/180.
PhotProb = []

def Double(m, params):
    return params[0]*(1.-tanh((m-params[1])/params[2]))*(1.-tanh((m-params[1])/params[3]))/4.

def Square(m, params):
    return (params[0]-params[1]*(m-21.0)**2)/(1.+exp((m-params[2])/params[3]))

def read_file_raw(fichier):
    """
    Reads in file 'fichier' and returns a list of lines.
    No processing is done on any line.
    """
    iter = 0
    while True:
        try:
	    f = file(fichier,'r')
	except IOError:
	    print 'Cannot open', fichier
	    iter += 1
	    if iter > 5:
	        print 'Something is wrong, quitting!'
                lines = None
		break
	    fichier = raw_input('Please re-enter the file name: ')
	else:
	    lines = f.readlines()
	    f.close()
            break
    return lines

def FixMag(fich):
    lines = read_file_raw(fich)
    fo = open('zzzzzzzz', 'w')
    if lines:
        for line in lines:
            if line[0] == '#':
                fo.write(line)
            else:
                words = line.split()
                if words[21] == 'g':
                    fo.write("%s  %5.2f %s %6.3f %7.3f %5.2f  g  %s\n" %(line[0:16], float(words[2])+0.0, line[24:163], float(words[18]), float(words[19]), float(words[20])+0.0, words[22]))
                if words[21] == 'r':
                    fo.write("%s  %5.2f %s %6.3f %7.3f %5.2f  g  %s\n" %(line[0:16], float(words[2])+0.7, line[24:163], float(words[18]), float(words[19]), float(words[20])+0.7, words[22]))
                if words[21] == 'R':
                    fo.write("%s  %5.2f %s %6.3f %7.3f %5.2f  g  %s\n" %(line[0:16], float(words[2])+0.8, line[24:163], float(words[18]), float(words[19]), float(words[20])+0.8, words[22]))
    fo.close()
    os.rename('zzzzzzzz', fich)

if __name__ == "__main__":
    FixMag('CFEPS-Classical.detections')
    FixMag('CFEPS-Plutino.detections')
    FixMag('CFEPS-DetachedOuter.detections')
    FixMag('CFEPS.detections')
    FixMag('CFEPS-Inner.detections')
    FixMag('CFEPS-noclassyet.detections')
    FixMag('CFEPS-Resonant.detections')
    FixMag('CFEPS-Scattering.detections')
    FixMag('CFEPS-SDO.detections')
