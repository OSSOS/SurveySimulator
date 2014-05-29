#! /usr/bin/env python
"""
Adds efficiency of detection at end of each line for a detected object.
"""

from numpy import tanh

def get_file_list(rep, pat):
    """Returns a list of files from directory 'rep', matching pattern 'pat'."""
    import re, os
    
    _sre_pat = re.compile(pat)
    _os_cwd = os.path.abspath(os.getcwd())
    os.chdir(rep)
    lsrep = os.listdir('.')
    lsrep = filter(_sre_pat.search, lsrep)
    lsrep.sort()
    lsrep = [os.getcwd()+'/'+e for e in lsrep]
    os.chdir(_os_cwd)
    
    return lsrep

def Efficiency(mag, eff):
    """
    Compute efficiency of detection for magnitude 'mag' from parameters in eff.
    """
    return eff[0]/4.*(1.-tanh((mag-eff[1])/eff[2]))*(1.-tanh((mag-eff[1])/eff[3]))

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

def read_eff(fichier):
    """
    Reads in file 'fichier' and record informations from efficiency file.
    """
    lines = read_file_raw(fichier)
    if (lines):
        for line in lines:
            if line[0:14] == 'double_param= ':
                vals = line[14:]
                try:
                    eff = [float(vals.split()[0]), float(vals.split()[1]), float(vals.split()[2]), float(vals.split()[3])]
                except ValueError:
                    pass
    return eff

def read_write_file(fichier, eff):
    """
    Reads in detection file and writes it out, appending efficiency of detection
    of each object.
    """
    lines = read_file_raw(fichier)
    if (lines):
        fo = open(fichier+'.eff', 'w')
        for line in lines:
            if line[0:1] == '#':
                fo.write('%s' % (line))
            else:
                key = line.split()[0][0:3]
                mag = float(line.split()[1])
                if key == 'L3h':
                    mag -= 0.7
                if key == 'K02':
                    mag -= 0.8
                fo.write('%s  %5.3f\n' % (line[:-1], Efficiency(mag, eff[key])))

        fo.close()
        return fichier+'.eff'

if __name__ == '__main__':
    eff = {}
    rep = '.'
    pattern = '.+-smooth\.eff'
    files = get_file_list(rep, pattern)
    for i in range(len(files)):
        fichier = files[i].split('/')[-1]
        key = fichier[0:3]
        if key == 'pre':
            key = 'K02'
        eff[key] = read_eff(fichier)

    pattern = '^CFEPS.+detections$'
    files = get_file_list(rep, pattern)
    for i in range(len(files)):
        fichier = files[i].split('/')[-1]
        fich = read_write_file(fichier, eff)
        print 'Wrote file: ' + fich

