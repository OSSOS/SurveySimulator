# GiMeObj

The GiMeObj subroutine is the critical routine which contains the model of the 
outer Solar System population(s) which are being exposed to the Survey
Simulator (and thus, ultimately, compared to real detections).  

When the driver code calls this routine, the driver expects that GiMeObj
will return one outer Solar System object defined by:

1. An osculating BARYCENTRIC ECLIPTIC J2000 orbit 
   	  The orbital element set (a, e, i, long_node,arg_peri, M, JD) must be
   	used, where the mean anomaly M's value is given at epoch JD.
 	  Note that the propogation forward or backward of the object's
	position in Detos1() is done as an unperturbed barycentric 
        2-body problem.

2. An absolute H magnitude and set of colours in the major filters.  
	  The H-mag can be specified in any filter (here called 'x') the
   	user wishes. The Bowell 'HG' formalism will be used to compute
	apparent magnitudes, so the 'G' parameter of Bowell is also supplied
  	as variable gb.
	  The user MUST fill the array with colours which at a minimum covers
	all filters used by the survey blocks that the user incorporates
  	into the survey they are simulating.  If all the detection blocks
	were in the 'g' filter, and the user specified H-magnitudes in the
        same filter (that is, x=g-band), then one could run a simulation
        where all colour terms could be zeroed.  As a next example, if the 
        detection imaging was all done in g or r band, then the user could
        specify H_g mags and provide only g-g=0 (color(1) below) and
        r-g = -0.7 (color(2) below).  One could even have 'x' be in a
        band which is none of the nine pre-defined bands (see below).
        Thus, the user only need specify the colors for the filters used
        in the detection surveys.
	  Note that the Survey Simulator Detos1() routine will RETURN a
  	magnitude in the band 'x' selected by the user.

3. A light curve amplitude, period, and phase 
	The user may turn this off by setting the amplitude to zero (the period
        must not be set to zero!).  The phase is interpreted as the rotational 
        phase at the same instant as the orbital element epoch (that is, the JD 
        at which M is given). The amplitude defined in the code is the 
        peak-to-peak amplitude.

NOTE: Fortran logical unit numbers 7 to 19 are reserved to use by `Driver.f` and
        `SurveySubs.f` and MUST NOT be used by the `GiMeObj` subroutine or any other 
    subroutine  that you may add to the driver. Please use logical unit numbers
	starting from 20.

----

## GiMeObj API

The list of arguments (arg_list_1 above) that the driver sends to GiMeObj is:

    (filena, seed, a, e, inc, node, peri, M, epoch, h, color, gb, ph,
     period, amp, comment, nchar, ierr)

with:

### INPUT
    filena: name of the file to be read in by GiMeObj the first time it 
            is called (CH)
    seed  : Random number generator seed (I4)

### OUTPUT
    a     	: semimajor axis (R8)
    e     	: eccentricity (R8)
    inc   	: Inclination with respect to J2000 ecliptic [rad] (R8)
    node  	: Longitude of ascending node [rad] (R8)
    peri  	: Argument of perihelion [rad] (R8)
    M     	: Mean anomaly [rad] (R8)
    epoch 	: epoch for M (and rotational phase below), in Julian Day (R8)
    h     	: absolute magnitude of object in band filter "x" (R8)
    color 	: array of colors "y-x", where the index of "y" is as
            	  described in detos1 (10*R8)
                   colors(1) : g-x
                   colors(2) : r-x
                   colors(3) : i-x
                   colors(4) : z-x
                   colors(5) : u-x
                   colors(6) : V-x
                   colors(7) : B-x
                   colors(8) : R-x
                   colors(9) : I-x
    gb    	: opposition surge factor G, Bowell formalism (R8)
    ph    	: phase of lightcurve at epoch [rad] (R8)
    period	: period of lightcurve [day] (R8)   CANNOT SET TO ZERO
    amp   	: peak-to-peak amplitude of lightcurve [mag] (R8)  
                  CAN  SET TO ZERO
    comment  	: user specified string containing whatever the user wants
		  (CHAR*100); can be empty.
    nchar	: number of characters in the comment string that should be
    		  printed out in output files if the object is detected;
                  maximum of 100 (I4)
    ierr  	: return code
                     0 : GiMeObj does not diagnose any errors, normal return
                         value
                   100 : end of model, exit after checking this object
                   -10 : could not get all orbital elements, skip object
                   -20 : something went grossly wrong, should quit

Normally the Driver itself terminates the simulator (based on enough tracked
detections having been gathered), but the return codes here allow the user the
option of having GiMeObj return codes to tell the driver to terminate.

---
## Other subroutines 

### Detos1		

Attempt to DETect 1 Outer Solarsystem object

The `Detos1` subroutine is the heart of the Survey Simulator.  

The driver calls `GiMeObj` and then asks `Detos1` : "Would this object be
seen by the surveys which are being simulated?"  More precisely, it says:
Given the detection efficiencies and pointing history of all blocks of
the survey, and taking into account the probabilistic nature of detection
(especially for faint objects), is this object in the field coverage and
detected by any block of the survey?  If so, where and when was it detected 
how bright, and was it tracked to a high precision orbit?

The list of arguments (arg_list_2 above) for Detos1 is

    (a, e, inc, node, peri, mt0, jday, hx, color, gb, ph, period, amp, surnam,
     seed, flag, ra, dec, d_ra, d_dec, r, delta, mi_int, m_rand, eff, isur, mt,
     jdayp, ic, surna, h_rand)

with:

#### INPUT (from GiMeObj)
    a     : Semi-major axis [AU] (R8)
    e     : Eccentricity (R8)
    inc   : Inclination [rad] (R8)
    node  : Longitude of node [rad] (R8)
    peri  : Argument of perihelion [rad] (R8)
    mt0   : Mean anomaly [rad] (R8)
    jday  : Reference time of the orbital elements [JD] (R8)
    hx    : Absolute magnitude of the object in the user's specified 'x' band
            (R8)
    color : Array of colors (10*R8)
               colors(1) : g-x
               colors(2) : r-x
               colors(3) : i-x
               colors(4) : z-x
               colors(5) : u-x
               colors(6) : V-x
               colors(7) : B-x
               colors(8) : R-x
               colors(9) : I-x
    gb    : opposition surge factor G, Bowell formalism (R8)
    ph    : phase of lightcurve at epoch jday [rad] (R8)
    period: period of lightcurve [day] (R8)
    amp   : peak-to-peak amplitude of lightcurve [mag] (R8)
    surnam: Survey directory name (CH10)

#### OUTPUT

    seed  : Random number generator seed (I4)
    flag  : Return flag (I4): 
		0: not found 
		1: found, but not tracked
            	2: found and tracked
    ra    : Right ascension at detection [rad] (R8)
    dec   : Declination at detection [rad] (R8)
    d_ra  : Right ascension rate [rad/day] (R8)
    d_dec : Declination rate [rad/day] (R8)
    r     : Sun-object distance [AU] (R8)
    delta : Earth-object distance [AU] (R8)
    m_int : Intrinsic apparent magnitude, in x-band (R8); This is computed from
            the absolute magnitude, the distance from the Sun and Earth and the
            phase angle (Bowel formalism); The returned value is given in the
            user defined 'x' filter
    m_rand: Averaged randomized magnitude, in x-band (R8); This is computed from
            the intrinsic apparent magnitude and a random uncertainty added
            according to Gaussian noise whose amplitude and center are
            determined from the <mag_error> parameters in the efficiency files;
            The returned value is given in the user defined 'x' filter
    eff   : Efficiency of detection of object (function of mag and survey) (R8)
    isur  : Identification number of the survey block the object was in (I4)
    mt    : Mean anomaly at discovery [rad] (R8)
    jdayp : Time of discovery [JD] (R8)
    ic    : Index of color used for survey (I4)
    surna : Detection survey name, for information (CH10) 
    h_rand: Absolute randomized magnitude (R8)

---

