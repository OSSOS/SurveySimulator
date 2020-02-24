import ssimTools
import numpy as np
import random
drad = np.pi/180.0


class Var:
    # Class to hold the initial variables

    H = None
    Hi = 0
    Hl = 0

    e = None
    ei = 0
    el = 0
    firste = True

    inc = None
    inci = 0
    incl = 0
    firstinc = True

    ai = 0
    a = 0
    al = 0
    firsta = True

    resamp = None
    resampi = 0
    resampl = 0
    firstresamp = True

    resamp2 = None
    resampi2 = 0
    resampl2 = 0
    firstresamp2 = True

    fdraw = True
    fsize = None
    rnd = None

    """
    def __init__(self, H = None, Hi = 0, Hl=0, 
        e = None, ei = 0, el = 0, firste = True, 
        inc = None, inci = 0, incl = 0, firstinc=True, 
        ai = 0, a = 0, al = 0, firsta = True, 
        resamp = None, resampi = 0, resampl = 0, firstresamp = True, 
        resamp2 = None, resampi2 = 0, resampl2 = 0, firstresamp2 = True, 
        fdraw = True, fsize = None, rnd = None):

        self.H = H
        self.Hi = Hi
        self.Hl= Hl

        self.e = e
        self.ei = ei
        self.el = el
        self.firste = firste

        self.inc = inc
        self.inci = inci
        self.incl = incl
        self.firstinc= firstinc

        self.ai = ai
        self.a = a
        self.al = al
        self.firsta = firsta

        self.resamp = resamp
        self.resampi = resampi
        self.resampl = resampl
        self.firstresamp = firstresamp

        self.resamp2 = resamp2
        self.resampi2 = resampi2
        self.resampl2 = resampl2
        self.firstresamp2 = firstresamp2

        self.fdraw = fdraw
        self.fsize = fsize
        self.rnd = rnd
    """


def setrand(seed):
    # Allow the seed to be set for repeatability
    Var.rnd = random.Random(seed)


def drawaold(exp):

    size = 10**6

    ind = Var.ai
    if Var.firsta is False and ind != Var.al:
        a = Var.a[ind]
        Var.ai += 1
        return a

    amin = 150
    amax = 1000

    if Var.firsta:
        Var.xpa = np.arange(amin, amax, 1)
        y = (Var.xpa ** exp).cumsum()
        Var.ya = y / max(y)

        Var.firsta = False

    rv2 = np.random.rand(size)
    Var.a = np.interp(rv2, Var.ya, Var.xpa)
    Var.a = Var.a[Var.a > amin]
    Var.al = len(Var.a)
    Var.ai = 1

    return Var.a[0]


def drawinc(mu, sigma):
    size = 10**6

    ind = Var.inci
    if Var.firstinc is False and ind != Var.incl:
        inc = Var.inc[ind]
        Var.inci += 1
        return inc

    incmin = 0
    incmax = 90*np.pi/180

    sigma = sigma*np.pi/180  # e.g. 6.9
    mu = mu*np.pi/180  # e.g. 19.3

    if Var.firstinc:
        Var.xpinc = np.arange(0, incmax, 0.01)
        y = (np.sin(Var.xpinc) * np.exp(-(Var.xpinc - mu) ** 2 / (2 * sigma ** 2))).cumsum()
        Var.yinc = y / max(y)

        Var.firstinc = False

    rv2 = np.random.rand(size)
    Var.inc = np.interp(rv2, Var.yinc, Var.xpinc)
    Var.inc = Var.inc[Var.inc > incmin]
    Var.incl = len(Var.inc)
    Var.inci = 0

    return Var.inc[0]


def drawecc(ec, ew):
    size = 10**6

    ind = Var.ei
    if Var.firste is False and ind != Var.el:
        e = Var.e[ind]
        Var.ei += 1
        return e

    emin = 0.0
    emax = 1.0

    if Var.firste:
        Var.xpe = np.arange(emin, emax, 0.01)
        y = (np.exp(-(Var.xpe - ec) ** 2 / (2 * ew ** 2))).cumsum()
        Var.ye = y / max(y)

        Var.firste = False

    rv2 = np.random.rand(size)
    Var.e = np.interp(rv2, Var.ye, Var.xpe)
    Var.e = Var.e[Var.e > emin]
    Var.el = len(Var.e)
    Var.ei = 0

    return Var.e[0]


def drawresamp(loamp, midamp, hiamp):
    size = 10**6

    ind = Var.resampi
    if Var.firstresamp is False and ind != Var.resampl:
        resamp = Var.resamp[ind]
        Var.resampi += 1
        return resamp

    if Var.firstresamp:
        Var.xpresamp = np.arange(loamp, hiamp, 0.01)
        y1 = (1. / (midamp-loamp) * Var.xpresamp[Var.xpresamp < midamp] - loamp / (midamp - loamp))
        y2 = (-1. / (hiamp-midamp) * Var.xpresamp[Var.xpresamp >= midamp] + hiamp / (hiamp - midamp))
        y = np.append(y1, y2).cumsum()
        Var.yresamp = y / max(y)

        Var.firstresamp = False

    rv2 = np.random.rand(size)
    Var.resamp = np.interp(rv2, Var.yresamp, Var.xpresamp)
    Var.resamp = Var.resamp[Var.resamp > loamp]
    Var.resampl = len(Var.resamp)
    Var.resampi = 0

    return Var.resamp[0]


def drawresamp2(loamp2, midamp2, hiamp2):
    size = 10**6

    ind = Var.resampi2
    if Var.firstresamp2 is False and ind != Var.resampl2:
        resamp2 = Var.resamp2[ind]
        Var.resampi2 += 1
        return resamp2

    if Var.firstresamp2:
        Var.xpresamp2 = np.arange(loamp2, hiamp2, 0.01)
        y1 = (1. / (midamp2-loamp2) * Var.xpresamp2[Var.xpresamp2 < midamp2] - loamp2 / (midamp2 - loamp2))
        y2 = (-1. / (hiamp2-midamp2) * Var.xpresamp2[Var.xpresamp2 >= midamp2] + hiamp2 / (hiamp2 - midamp2))
        y = np.append(y1, y2).cumsum()
        Var.yresamp2 = y / max(y)

        Var.firstresamp2 = False

    rv2 = np.random.rand(size)
    Var.resamp2 = np.interp(rv2, Var.yresamp2, Var.xpresamp2)
    Var.resamp2 = Var.resamp2[Var.resamp2 > loamp2]
    Var.resampl2 = len(Var.resamp2)
    Var.resampi2 = 0

    return Var.resamp2[0]


def gimeobj(fname=None):
    """
    GiMiObj returns an object with orbital elements and observation characteristics.
    GiMiObj takes in a file of orbital elements (the KRQ q200 file in this case),
    selects an object, "fuzzes" its orbital parameters
    assigns an H magnitude and a colour and then returns it.
    gimeobj is currently configured to choose n:1 resonators,
    assign an H magnitude and a colour, and then return all of that.

    Inputs:
    fname - a file in this directory

    Output:
    a - semimajor axis
    e - eccentricity
    inc - inclination (rad)
    node - longitude of the ascending node (rad)
    peri - argument of pericentre (rad)
    M - mean anomaly (rad)
    h - H magnitude in field of survey
    color - array of colour conversion values (see below)
    gb - opposition surge effect
    ph - phase angle
    period - light curve period
    amp - light curve amplitude
    """

    if Var.fdraw:
        # Read in H distribution parameters plus one other deprecated value ....
        Var.alpha_faint, Var.alpha, Var.hmax, Var.hmin, Var.hbreak, Var.contrast = \
            np.genfromtxt('H_dist_vars.txt', unpack=True)

        # No need to reload the file
        Var.fdraw = False

    rg = Var.rnd

    # ------------- Determine orbital parameters -----------------------

    def draw_candidate(ares, loamp, midamp, hiamp):
        # choose a within 0.5 of 30.1*(2/1)^(2/3) = 47.8
        a_draw = rg.random() * 1.0 + ares

        # choose e from gaussian with center ec and width ew
        e_draw = 1.-30.1/ares

        # choose inc from sin(i)*gaussian with center mu and width sigma
        inc_draw = 0.001

        # choose libration amplitude, then phi sinusoidally from that lib amp
        resamp_draw = drawresamp(loamp, midamp, hiamp)

        return a_draw, e_draw, inc_draw, resamp_draw

    def draw_candidate2(ares, loamp, midamp, hiamp):
        # choose a within 0.5 of 30.1*(2/1)^(2/3) = 47.8
        a_draw = rg.random() * 1.0 + ares

        # choose e from gaussian with center ec and width ew
        e_draw = 1.-30.1/ares

        # choose inc from sin(i)*gaussian with center mu and width sigma
        inc_draw = 0.001

        # choose libration amplitude, then phi sinusoidally from that lib amp
        resamp_draw = drawresamp2(loamp, midamp, hiamp)

        return a_draw, e_draw, inc_draw, resamp_draw

    # different for sym or asym
    if rg.random() < 0.3:  # symmetric
        a, e, inc, resamp = draw_candidate(47.8, 0, 5, 10)
        # choose value of phi sinusoidally from resamp
        phi = np.pi + np.sin(2.*np.pi*rg.random())*resamp/180.*np.pi
    else:  # asymmetric
        a, e, inc, resamp = draw_candidate2(47.8, 0, 5, 10)
        phi = 80./180.*np.pi + np.sin(2.*np.pi*rg.random())*resamp/180.*np.pi
        if rg.random() > 0.5:  # half in leading island, half in trailing island
            phi = 2.*np.pi-phi

    # Assign two angles randomly
    node = rg.random()*2*np.pi
    m = rg.random()*2*np.pi

    # choose value of phi sinusoidally from resamp (chosen above)
    # phi=np.pi +np.sin(2.*np.pi*rg.random())*resamp/180.*np.pi

    # assign last angle so that phi74 works
    lambda_n = 333.5078*np.pi/180.  # Neptune's mean longitude on 1 Jan 2013
    peri = (1./1.*(phi - 2.*m) - node + lambda_n) % (2.*np.pi)

    # Draw an H magnitude from a SPL, knee, or divot distribution
    # Uses class values as read in from file (see above)
    h = ssimTools.draw_h(Var.alpha, Var.hmax, Var.alpha_faint, Var.contrast, Var.hbreak, Var.hmin)

    epoch = 2456293.5
    gb = -0.12      # Opposition surge effect
    ph = 0.00       # Initial phase of lightcurve
    period = 0.60   # Period of lightcurve
    amp = 0.00      # Amplitude of lightcurve (peak-to-peak)
    #     color : Array of colors (10*R8)
    #                colors[0] : g-x
    #                colors[1] : r-x
    #                colors[2] : i-x
    #                colors[3] : z-x
    #                colors[4] : u-x
    #                colors[5] : V-x
    #                colors[6] : B-x
    #                colors[7] : R-x
    #                colors[8] : I-x  

    gmr = rg.gauss(0.7, 0.2)  # g-r colour sampled from gaussian with mu 0.7 and std 0.2 (from CFEPS data)
    color = [gmr, 0.0, 0.0, 0.0, 0.0, 0.0,  0.0, -0.1, 0.0, 0.0]
    # This color array is set-up for r band (note r-x = 0.0). 
    # There is a required 10th element which currently does not correspond to a colour.
    # It is there for future expansion, or so I'm told.

    return a, e, inc, node, peri, m, epoch, h, color, gb, ph, period, amp, resamp

    # The call to SurveySubs.detos1 for reference
    # SurveySubs.detos1(a,e,inc,node,peri,M,epoch,h,color,gb,ph,period,amp,survey_dir,seed)
