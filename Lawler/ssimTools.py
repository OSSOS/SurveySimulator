import numpy as np
from random import gauss
import random
drad = np.pi/180.0

class iv:
# Class to hold the initial variables (iv)
    y = None
    xp = None
    first = True
    l1 = None

    H = None
    Hi = 0
    Hl = 0
    rnd = None

def setRand(seed):
    iv.rnd = random.Random(seed)
    np.random.seed(seed) 
   
def drawH(alpha, hmax, alpha_faint=None, contrast=1, hbreak=None, hmin=1):

    """Compute and assign and H-magnitude from a so-called single

    power-law, knee, or divot H-magnitude distribution.

    When provided a slope alpha and a faint-side maximum H-magnitude

    (hmax), a H-magnitude is drawn randomly from the distribution

                       dN/dH propto 10**(alpha H)

    in the range hmin = 1 to hmax. Specify hmin to change the bright-end.

    Specifying an hbreak and alpha_faint will draw from a knee distribution

    Specifying an hbreak, alpha_faint and contrast will draw from a divot

    distrubtion as in Shankman et al. 2013

    e.g.

    ---Single Power Law---

    drawH(0.8,13)

    will draw an H-magnitude from the appropriate distribution such that

    H [1,13]

    drawH(0.8,13,hmin=5) 

    will draw an H-magnitude such that H [5,13]

    ---Knee---

    To draw from a knee distribution specify hbreak and alpha_faint

    drawH(0.8, 13, hbreak=9, alpha_faint = 0.5)

    This will draw an H-magnitude from a distrubtion that breaks at H=9

    from a slope of 0.8 to a slope of 0.5. hmin can also be specified here.


    ---Divot---

    To draw from a divot (see Shankman et al 2013), specify hbreak,

    alpha_faint, and the contrast value. Contrasts should be > 1.

    hmin can also be specified.

    drawH(0.8, 13, hbreak=9, alpha_faint = 0.5, contrast = 23)

    """
    # This approach generates the distribution as an x,y pair. For a SPL this is trivially a straight line.
    # For knees and divots, the pre and post break parts must be made separately then properly joined at the
    # break in a differental way and then turned into a cumulative distribution. 
    #
    # Once this is done, a large number of H objects are drawn from the distribution using interpolation.
    # The x and y axis are switched, a set of random numbers between 0 and 1 (for y) is generated and then for each
    # random number, an H value (x) is interpolated. These values are then stored in an array in iv. Subsequent calls
    # to the program will just return the next H already drawn. If the number of calls exceeds the number initially drawn,
    # a new large sample of H values will be drawn. This is done for speed.


    size = 10**6 # The size of the H array that will be generated

    i = iv.Hi # The current index location in the H array
    # If this is not the first draw and if the index is not at the limit then return the next H value in the array
    # and increase the counter
    if iv.first == False and i!= iv.Hl:
        H = iv.H[i]
        iv.Hi += 1
        return H
    
    # Simple error handling for silly input
    if hmax < hbreak:
        raise ValueError('hmax must be greater than hbreak')
    if hmax < hmin:
        raise ValueError('hmax must be greater than hmin')
    if hbreak is not None and hbreak < hmin:
        raise ValueError('hbreak must be greater than hmin')

    # Avoid singularity for alpha = 0
    alpha = 0.0000000001 if alpha == 0 else alpha

    # Set alpha_faint to alpha for the case of a single power-law
    alpha_faint = alpha if alpha_faint is None else alpha_faint

    # Avoid singularity for alpha_faint = 0
    alpha_faint = 0.0000000001 if alpha_faint == 0 else alpha_faint

    # Set hbreak to be the maximum H for the case of a single power-law
    hbreak = hmax if hbreak is None else hbreak


    # This part constructs the cumulative size distribution.
    # **** Note that the size distribution cannot be changed once set. ****
    if iv.first==True:
        
        # Generate the bright side of the size distribution
        xpb = np.arange(hmin-(hmax-hmin)*0.01,hbreak,0.001)
        yb =10**(alpha*xpb)*alpha*np.log10(10)

        # Generate the faint side of the size distribution
        xpf = np.arange(hbreak,hmax,0.001)
        yf = 10**(alpha_faint*xpf)

        # Merge the x axes into one array 
        xp = np.concatenate([xpb, xpf])
        # Pad first element with the minimum H value (essentially setting it up properly with correct bin edges)
        iv.xp = np.zeros(len(xp)+1)+hmin
        iv.xp[1:] = xp
        

        # Coefficient to merge the faint and bright y arrays
        yf = yf*( alpha/alpha_faint * 10**(alpha*hbreak)/(10**(alpha_faint*hbreak))/contrast)*alpha_faint*np.log10(10)
        
        # Merge them and then cumulate and normalzie
        y = np.concatenate([yb,yf]).cumsum()
        y = y/max(y)
        
        #Pad cumulative distribution with a 0 as all cumulatives should start there
        iv.y = np.zeros(len(y) +1)
        iv.y[1:] = y

        # This only needs to be done once. Set first to False
        iv.first = False


    # Generate a set of random numbers between 0 and 1 for the interpolation
    rv2 = np.random.rand(size)
    # Flip x and y and interpolate it for input y's, thus given an array of x's (H values)
    iv.H = np.interp(rv2, iv.y, iv.xp)
    # Store the initial length of the H array
    iv.Hl = len(iv.H)
    # The counter is set to 1 as the first element will already be returned
    iv.Hi = 1
    
    # The first element is returned
    return iv.H[0]


def fuzz(var, pct, type=None):
    # Fuzzes a variable var randomly by up +- percent pct as a percent
    # or if a type is specified (anything provided that is not "None"  will do) then it fuzzes by up to +- the value
    # e.g. a = fuzz(a, 0.1) fuzzes a by up to +- 10 pct
    #      a = fuzz(a, 10) also fuzzes a by up to +- 10 pct
    # but  inc = fuzz(inc, 1, type = 'abs') will fuzz i by up to +- 1
    rnd = iv.rnd
    pct = pct/100.0 if (pct > 1.0 and type is None) else pct
    var = (var*(1.0 + pct*(2.0*rnd.random()-1.0))) if type is None else (var + (2.0*rnd.random()-1.0)*pct)
    return var



# ***** Handle all outputs in the OSSOS format *****

def detFile(fname, seed):
    # Open detection file and write header
    f_detect = open(fname, 'w')

    f_detect.write ('# Seed: %10d\n#\n' % (seed))
    f_detect.write ('#   a      e        i        node     peri     Manom      resamp       q        r        M      m_rand H_rand color flag delta    m_int   H_int eff   RA(H)     DEC    delta_ra delt_dec Surv.  Comments\n#\n')
    f_detect.write ('# flag: >0: detected; >2: characterized; 0 mod(2): tracked\n')
    f_detect.write ('# Survey: name of the block\n')
    f_detect.write ('# delta_ra: distance from center of pointing [arcsec]\n')
    f_detect.write ('# delt_dec: distance from center of pointing [arcsec]\n#\n')


def trackFile(fname):
    # Open tracked detection file (no header)
    f_track = open(fname, 'w')
    f_track.write ('#   a      e        i       node      peri     Manom      resamp      q        r        M      m_rand H_rand color Comment\n#\n')
    f_track.close()

def detWrite(fname, a, e, inc, node, peri, M, resamp, r, mt, m_rand, h_rand, color, ic, flag, delta, m_int, h, eff, ra, dec, d_ra, d_dec, surna, comments):
    f_detect = open(fname, 'a')
    f_detect.write('%8.3f %6.3f %8.3f  %8.1f  %8.1f  %8.1f %8.1f %8.3f %8.3f %8.3f %8.3f %6.2f %5.2f %2d %8.3f %8.3f %6.2f %4.2f %8.5f %8.4f %8.5f %8.5f %6s %s\n' % (a, e, inc/drad, node/drad, peri/drad, M/drad, resamp, a*(1.-e), r, mt/drad, m_rand, h_rand, color[ic-1], flag, delta, m_int, h, eff, ra/drad/15., dec/drad, d_ra/drad*3600./24., d_dec/drad*3600./24., surna, comments))
    f_detect.close()

def trackWrite(fname, a, e, inc, node, peri, M, resamp, r, mt, m_rand, h_rand, color, ic,  comments):
    f_track = open(fname, 'a')
    f_track.write("%8.3f %6.3f %8.3f %8.1f  %8.1f  %8.1f %8.1f  %8.3f %8.3f %8.3f %8.3f %6.2f %5.2f %s\n" % (a, e, inc/drad, node/drad, peri/drad, M/drad, resamp, a*(1.-e), r, mt/drad, m_rand, h_rand, color[ic-1], comments))
    f_track.close()

def detSuffix(fname, n_iter, n_hits, n_track):
    f_detect = open(fname,'a')
    f_detect.write ('# Total number of objects:   %11d\n' % (n_iter))
    f_detect.write ('# Number of detections:      %7d\n' % (n_hits))
    f_detect.write ('# Number of tracked objects: %7d\n' % (n_track))
    f_detect.close()
