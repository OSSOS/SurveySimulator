import random

import numpy as np

drad = np.pi/180.0


class Vars:
    # Class to hold the initial variables (iv)
    y = None
    xp = None
    first = True
    l1 = None

    H = None
    Hi = 0
    Hl = 0
    rnd = None


def setrand(seed):
    Vars.rnd = random.Random(seed)
    np.random.seed(seed)


def draw_h(alpha, hmax, alpha_faint=None, contrast=1, hbreak=None, hmin=1):

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

    draw_h(0.8,13)

    will draw an H-magnitude from the appropriate distribution such that

    H [1,13]

    draw_h(0.8,13,hmin=5)

    will draw an H-magnitude such that H [5,13]

    ---Knee---

    To draw from a knee distribution specify hbreak and alpha_faint

    draw_h(0.8, 13, hbreak=9, alpha_faint = 0.5)

    This will draw an H-magnitude from a distrubtion that breaks at H=9

    from a slope of 0.8 to a slope of 0.5. hmin can also be specified here.


    ---Divot---

    To draw from a divot (see Shankman et al 2013), specify hbreak,

    alpha_faint, and the contrast value. Contrasts should be > 1.

    hmin can also be specified.

    draw_h(0.8, 13, hbreak=9, alpha_faint = 0.5, contrast = 23)

    """
    # This approach generates the distribution as an x,y pair. For a SPL this is trivially a straight line.
    # For knees and divots, the pre and post break parts must be made separately then properly joined at the
    # break in a differential way and then turned into a cumulative distribution.
    #
    # Once this is done, a large number of H objects are drawn from the distribution using interpolation. The x and y
    # axis are switched, a set of random numbers between 0 and 1 (for y) is generated and then for each random
    # number, an H value (x) is interpolated. These values are then stored in an array in iv. Subsequent calls to the
    # program will just return the next H already drawn. If the number of calls exceeds the number initially drawn,
    # a new large sample of H values will be drawn. This is done for speed.

    size = 10**6  # The size of the H array that will be generated

    i = Vars.Hi  # The current index location in the H array
    # If this is not the first draw and if the index is not at the limit then return the next H value in the array
    # and increase the counter
    if Vars.first is False and i != Vars.Hl:
        h = Vars.H[i]
        Vars.Hi += 1
        return h

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
    if Vars.first:

        # Generate the bright side of the size distribution
        xpb = np.arange(hmin-(hmax-hmin)*0.01, hbreak, 0.001)
        yb = 10**(alpha*xpb)*alpha*np.log(10)

        # Generate the faint side of the size distribution
        xpf = np.arange(hbreak, hmax, 0.001)
        yf = 10**(alpha_faint*xpf)

        # Merge the x axes into one array 
        xp = np.concatenate([xpb, xpf])
        # Pad first element with the minimum H value (essentially setting it up properly with correct bin edges)
        Vars.xp = np.zeros(len(xp) + 1) + hmin
        Vars.xp[1:] = xp

        # Coefficient to merge the faint and bright y arrays
        yf = yf*(alpha/alpha_faint * 10**(alpha*hbreak)/(10**(alpha_faint*hbreak))/contrast)*alpha_faint*np.log(10)

        # Merge them and then cumulate and normalzie
        y = np.concatenate([yb, yf]).cumsum()
        y = y/max(y)

        # Pad cumulative distribution with a 0 as all cumulatives should start there
        Vars.y = np.zeros(len(y) + 1)
        Vars.y[1:] = y

        # This only needs to be done once. Set first to False
        Vars.first = False

    # Generate a set of random numbers between 0 and 1 for the interpolation
    rv2 = np.random.rand(size)
    # Flip x and y and interpolate it for input y's, thus given an array of x's (H values)
    Vars.H = np.interp(rv2, Vars.y, Vars.xp)

    # Store the initial length of the H array
    Vars.Hl = len(Vars.H)
    # The counter is set to 1 as the first element will already be returned
    Vars.Hi = 1

    # The first element is returned
    return Vars.H[0]

# ***** Handle all outputs in the OSSOS format *****


def detfile(fname, seed):
    # Open detection file and write header
    f_detect = open(fname, 'w')

    f_detect.write(f'# Seed: {seed:10d}\n#\n')
    f_detect.write('# flag: >0: detected; >2: characterized; 0 mod(2): tracked\n')
    f_detect.write('# Survey: name of the block\n')
    f_detect.write('# delta_ra: distance from center of pointing [arcsec]\n')
    f_detect.write('# delt_dec: distance from center of pointing [arcsec]\n#\n')
    f_detect.write(f'#{"a [AU]":>7s} {"e":>6s} {"i [°]":>8s} {"Ω [°]":>8s} {"ω [°]":>8s} {"M [°]":>8s} '
                   f'{"ResAmp [°]":>10s} {"q [°]":>8s} {"Dist_* [AU]":>12s} {"M(t) [°]":>8s} {"MagR":>8s} '
                   f'{"H MagR":>6s} {"Color":>5s} {"flag":>5s} {"Dist_E [AU]":>12s} {"Mag_Intr":>8s} '
                   f'{"H_Intr":>6s} {"eff":>4s} {"RA [H]":>8s} {"Dec [°]":>8s} \u202F{"𝛿RA [″/hr]":>11s} '
                   f'\u200A{"𝛿DEC [″/hr]":>11s} {"Survey":>10s} {"Comments":>10s}\n#\n')
    # Special Unicode spaces (U+202F Narrow No-Break Space, and U+2009 Thin Space) are needed to ensure proper
    # spacing in the text file. This is because the double prime symbol used to denote arcseconds is thinner than
    # standard fixed width unicode characters.


def trackfile(fname):
    # Open tracked detection file and write header
    f_track = open(fname, 'w')
    f_track.write(f'#{"a [AU]":>7s} {"e":>6s} {"i [°]":>8s} {"Ω [°]":>8s} {"ω [°]":>8s} {"M [°]":>8s} '
                  f'{"ResAmp [°]":>10s} {"q [°]":>8s} {"Dist_* [AU]":>12s} {"M(t) [°]":>8s} {"MagR":>8s} '
                  f'{"H MagR":>6s} {"Color":>5s} {"Comments":>10s}\n#\n')
    f_track.close()


def detwrite(fname, a, e, inc, node, peri, m, resamp, r, mt, m_rand, h_rand, color, ic, flag,
             delta, m_int, h, eff, ra, dec, d_ra, d_dec, surna, comments):
    f_detect = open(fname, 'a')
    f_detect.write(f'{a:8.3f} {e:6.3f} {inc / drad:8.3f} {node / drad:8.1f} {peri / drad:8.1f} {m / drad:8.1f} '
                   f'{resamp:10.5f} {a * (1 - e):8.3f} {r:12.5f} {mt / drad:8.3f} {m_rand:8.3f} '
                   f'{h_rand:6.2f} {color[ic - 1]:5.2f} {flag:>5d} {delta:12.5f} {m_int:8.3f} '
                   f'{h:6.2f} {eff:4.2f} {ra / drad / 15:8.5f} {dec / drad:8.4f} {d_ra / drad * 3600 / 24:11.6f} '
                   f'{d_dec / drad * 3600 / 24:11.6f} {surna.decode("utf-8"):>10s} {comments:>10s}\n')
    f_detect.close()


def trackwrite(fname, a, e, inc, node, peri, m, resamp, r, mt, m_rand, h_rand, color, ic, comments):
    f_track = open(fname, 'a')
    f_track.write(f'{a:8.3f} {e:6.3f} {inc / drad:8.3f} {node / drad:8.1f} {peri / drad:8.1f} {m / drad:8.1f} '
                  f'{resamp:10.5f} {a * (1 - e):8.3f} {r:12.5f} {mt / drad:8.3f} {m_rand:8.3f} '
                  f'{h_rand:6.2f} {color[ic - 1]:5.2f} {comments:>10s}\n')
    f_track.close()


def detsuffix(fname, n_iter, n_hits, n_track):
    f_detect = open(fname, 'a')
    f_detect.write('#\n')
    f_detect.write(f'# Total number of objects:  {n_iter:11d}\n')
    f_detect.write(f'# Number of detections:     {n_hits:11d}\n')
    f_detect.write(f'# Number of tracked objects:{n_track:11d}\n')
    f_detect.close()
