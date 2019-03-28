from ssimTools import drawH, fuzz
import numpy as np
import random
drad = np.pi/180.0

class iv:
# Class to hold the initial variables (iv)
    H = None
    Hi = 0
    Hl=0
#    a = None

    e = None
    ei = 0
    el = 0
    firste = True

    inc = None
    inci = 0
    incl = 0
    firstinc=True

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

def setRand(seed):
    # Allow the seed to be set for repeatability
    iv.rnd = random.Random(seed)



def drawaold(exp):
    
    size = 10**6

    ind = iv.ai
    if iv.firsta == False and ind!= iv.al:
        a = iv.a[ind]
        iv.ai += 1
        return a

    amin = 150
    amax = 1000

    if iv.firsta==True:
       iv.xpa = np.arange(amin,amax, 1)
       y = (iv.xpa**(exp)).cumsum()
       iv.ya = y/max(y)

       iv.firsta = False

    rv2 = np.random.rand(size)
    iv.a = np.interp(rv2, iv.ya, iv.xpa)
    iv.a = iv.a[iv.a>amin]
    iv.al = len(iv.a)
    iv.ai = 1
    
    return iv.a[0]



def drawinc(mu,sigma):
    size = 10**6

    ind = iv.inci
    if iv.firstinc == False and ind!= iv.incl:
        inc = iv.inc[ind]
        iv.inci += 1
        return inc

    incmin = 0
    incmax = 90*np.pi/180

    sigma = sigma*np.pi/180 # e.g. 6.9     
    mu = mu*np.pi/180 #e.g. 19.3 

    if iv.firstinc==True:
       iv.xpinc = np.arange(0,incmax, 0.01)
       y = (np.sin(iv.xpinc)*np.exp(-(iv.xpinc-mu)**2/(2*sigma**2))).cumsum()
       iv.yinc = y/max(y)

       iv.firstinc = False

    rv2 = np.random.rand(size)
    iv.inc = np.interp(rv2, iv.yinc, iv.xpinc)
    iv.inc = iv.inc[iv.inc>incmin]
    iv.incl = len(iv.inc)
    iv.inci = 0

    return iv.inc[0]

def drawecc(ec,ew):
    size = 10**6

    ind = iv.ei
    if iv.firste == False and ind!= iv.el:
        e = iv.e[ind]
        iv.ei += 1
        return e

    emin = 0.0
    emax = 1.0
                                                                                                                                                                           
    if iv.firste==True:
       iv.xpe = np.arange(emin,emax, 0.01)
       y = (np.exp(-(iv.xpe-ec)**2/(2*ew**2))).cumsum()
       iv.ye = y/max(y)

       iv.firste = False

    rv2 = np.random.rand(size)
    iv.e = np.interp(rv2, iv.ye, iv.xpe)
    iv.e = iv.e[iv.e>emin]
    iv.el = len(iv.e)
    iv.ei = 0

    return iv.e[0]

def drawresamp(loamp,midamp,hiamp):
    size = 10**6

    ind = iv.resampi
    if iv.firstresamp == False and ind!= iv.resampl:
        resamp = iv.resamp[ind]
        iv.resampi += 1
        return resamp

                                                                                                                                                                          
    if iv.firstresamp==True:
       iv.xpresamp = np.arange(loamp,hiamp, 0.01)
       y1 = (1./(midamp-loamp)*iv.xpresamp[iv.xpresamp<midamp]-loamp/(midamp-loamp))
       y2 = (-1./(hiamp-midamp)*iv.xpresamp[iv.xpresamp>=midamp]+hiamp/(hiamp-midamp))
       y=np.append(y1,y2).cumsum()
       iv.yresamp = y/max(y)

       iv.firstresamp = False

    rv2 = np.random.rand(size)
    iv.resamp = np.interp(rv2, iv.yresamp, iv.xpresamp)
    iv.resamp = iv.resamp[iv.resamp>loamp]
    iv.resampl = len(iv.resamp)
    iv.resampi = 0

    return iv.resamp[0]

def drawresamp2(loamp2,midamp2,hiamp2):
    size = 10**6

    ind = iv.resampi2
    if iv.firstresamp2 == False and ind!= iv.resampl2:
        resamp2 = iv.resamp2[ind]
        iv.resampi2 += 1
        return resamp2

                                                                                                                                                                          
    if iv.firstresamp2==True:
       iv.xpresamp2 = np.arange(loamp2,hiamp2, 0.01)
       y1 = (1./(midamp2-loamp2)*iv.xpresamp2[iv.xpresamp2<midamp2]-loamp2/(midamp2-loamp2))
       y2 = (-1./(hiamp2-midamp2)*iv.xpresamp2[iv.xpresamp2>=midamp2]+hiamp2/(hiamp2-midamp2))
       y=np.append(y1,y2).cumsum()
       iv.yresamp2 = y/max(y)

       iv.firstresamp2 = False

    rv2 = np.random.rand(size)
    iv.resamp2 = np.interp(rv2, iv.yresamp2, iv.xpresamp2)
    iv.resamp2 = iv.resamp2[iv.resamp2>loamp2]
    iv.resampl2 = len(iv.resamp2)
    iv.resampi2 = 0

    return iv.resamp2[0]


def GiMeObj(fname):
    """
    GiMeObj is currently configured to choose n:1 resonators,
    assign an H magnitude and a colour, and then return all of that.

    Inputs:
    fname - a file in this director, need to mess with this

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



    if iv.fdraw==True:
        # Read in H distribution parameters plus one other deprecated value ....
        iv.alpha_faint, iv.alpha, iv.hmax, iv.hmin, iv.hbreak, iv.contrast = np.genfromtxt('H_dist_vars.txt', unpack=True)

        
        # No need to reload the file
        iv.fdraw = False

    
    rg = iv.rnd

    # ------------- Determine oribal parameters -----------------------

    def draw_candidate(ares,loamp,midamp,hiamp):
        #choose a within 0.5 of 30.1*(2/1)^(2/3) = 47.8
        a = rg.random()*(1.0) +ares
	
        #choose e from gaussian with center ec and width ew
        e = 1.-30.1/ares

        #choose inc from sin(i)*gaussian with center mu and width sigma
        inc = 0.001

        #choose libration amplitude, then phi sinusoidally from that lib amp
        resamp=drawresamp(loamp,midamp,hiamp)

        return a,e,inc,resamp

    def draw_candidate2(ares,loamp,midamp,hiamp):
        #choose a within 0.5 of 30.1*(2/1)^(2/3) = 47.8
        a = rg.random()*(1.0) +ares
	
        #choose e from gaussian with center ec and width ew
        e = 1.-30.1/ares

        #choose inc from sin(i)*gaussian with center mu and width sigma
        inc = 0.001

        #choose libration amplitude, then phi sinusoidally from that lib amp
        resamp=drawresamp2(loamp,midamp,hiamp)

        return a,e,inc,resamp

    #different for sym or asym
    if rg.random()<0.3:  #symmetric
      a, e, inc, resamp = draw_candidate(47.8,0,5,10)
      #choose value of phi sinusoidally from resamp
      phi=np.pi +np.sin(2.*np.pi*rg.random())*resamp/180.*np.pi
    else:  #asymmetric
      a, e, inc, resamp = draw_candidate2(47.8,0,5,10)
      phi=80./180.*np.pi +np.sin(2.*np.pi*rg.random())*resamp/180.*np.pi
      if rg.random()>0.5:  #half in leading island, half in trailing island
        phi=2.*np.pi-phi

    # Assign two angles randomly
    node = rg.random()*2*np.pi
    M = rg.random()*2*np.pi

    #choose value of phi sinusoidally from resamp (chosen above)
    #phi=np.pi +np.sin(2.*np.pi*rg.random())*resamp/180.*np.pi

    #assign last angle so that phi74 works
    lambdaN=333.5078*np.pi/180.  #Neptune's mean longitude on 1 Jan 2013
    peri = (1./1.*(phi - 2.*M) - node + lambdaN)%(2.*np.pi)

    # Draw an H magnitude from a SPL, knee, or divot distribution
    # Uses class values as read in from file (see above)
    h = drawH(iv.alpha, iv.hmax, iv.alpha_faint, iv.contrast, iv.hbreak, iv.hmin)


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

    gmr = rg.gauss(0.7, 0.2) # g-r colour sampled from gaussian with mu 0.7 and std 0.2 (from CFEPS data)
    color=[gmr, 0.0, 0.0, 0.0, 0.0, 0.0,  0.0, -0.1, 0.0, 0.0]
    # This color array is set-up for r band (note r-x = 0.0). 
    # There is a required 10th element which currently does not correspond to a colour.
    # It is there for future expansion, or so I'm told.

    return a, e, inc, node, peri, M, epoch, h, color,gb,ph,period,amp, resamp

    # The call to SurveySubs.detos1 for reference
    # SurveySubs.detos1(a,e,inc,node,peri,M,epoch,h,color,gb,ph,period,amp,survey_dir,seed)

