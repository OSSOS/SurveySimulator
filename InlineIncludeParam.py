#! /usr/bin/env python
"""
Inline the "include 'param.inc'" statements in SurveySubs.f.
"""

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
            try:
                fg = gzip.open(fichier+'.gz','r')
            except IOError:
                print 'Cannot open', fichier, ' nor ', fichier+'.gz'
                iter += 1
                if iter > 5:
                    print 'Something is wrong, quitting!'
                    lines = None
                    break
                fichier = raw_input('Please re-enter the file name: ')
            else:
                lines = fg.readlines()
                fg.close()
                break
        else:
	    lines = f.readlines()
	    f.close()
            break
    return lines

if __name__ == "__main__":
    fin = 'zzzz1'
    fout = 'zzzz0'
    lines = read_file_raw(fin)
    pas = read_file_raw('param.inc')
    if lines:
        fo = open(fout, 'w')
        for line in lines:
            cp = True
            ws = line.split()
            if (len(ws) == 2):
                if (((ws[0] == "include") or (ws[0] == "INCLUDE")) and (ws[1] == "'param.inc'")):
                    cp = False
                    for pa in pas:
                        fo.write(pa)
            if cp:
                fo.write(line)
        fo.close()
