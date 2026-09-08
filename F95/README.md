## Survey Simulator F95 based source

### Requirements

- gfortran
- make

### Install

e.g. `make Driver GIMEOBJ=InnerHotModel` 

links `InnerHotModel.f95` to `GeMeObj.f95` and then compiles `Driver`.

The `Makefile` builds `Driver` but first links a particular `GeMeObj` module to `GeMeObj.f95`. 

Users build `GeMeObj` subroutines and link them when building `Driver` then the simulator
will use that subroutine to get orbital elements and H magnitude for use in the simulation.

If you write MyGeMeObj.f95 model the build Driver as (leave off the .f95):

`make Driver GIMEOBJ=MyGeMeObj`

README.md (this file)  Describes "GiMeObj" routine
- InnerHotModel.f95 : Example "GiMeObj" routine to generate objects via parametric  prescription of inner belt
- ReadModelFromFile.f95 : Example "GiMeObj" routine to select objects from a model file. See Model File Format section for details.
- Makefile : Makefile to build executable
- \example : Contains files to run an example and test code


### COMPILING
  One can compile the program directly. Because Driver.f 'includes' a
  file 'GiMeObj.f' containing the model definition, one must first create 
  a symbolic link pointing to the actual file "InnerHotModel.f":

    ln -s InnerHotModel.f9 GiMeObj.f9
    gfortran -O3 -o SurveySimulator Driver.f9

or to read a model from an input file:

    ln -s ReadModelFromFile.f9 GiMeObj.f9
    gfortran -O3 -o SurveySimulator Driver.f9


### RUNNING

  `Driver` reads six parameters from the standard input:

  - the seed for the random number generator (integer)
  - a number <n> to control how long we run:
      <n> > 0: maximum number of simulated tracked detections
      <n> < 0: -maximum number of iterations (i.e. the number of calls to
               GiMeObj)
      <n> = 0: run until the 'model' decides to stop
  - the name of the directory containing the characterized survey blocks
    (see [`../docs/`](../docs/) for characterization formats and a mini example;
    full survey packages are distributed in the sibling
    [`SurveySimulator-Data`](../../SurveySimulator-Data/) tree;
    Driver test fixtures are under `tests/Surveys/`)
  - the name of the model input file (a file GiMeObj will read in)
  - the name of the output file where the detected objects will be listed
  - the name of the output file where the detected and tracked objects will be
    listed


  Normally one creates a file (e.g. `Driver.in`) with one line for each value, 
  and executes `Driver` like so:

    Driver < Driver.in

  An example of `Driver.in` can be found in the example subdirectory.

  Execution generates two output files with the user specified names:
  - The first (`SimulDetect.dat` in the example `Driver.in`) contains all
    detected objects. The meaning of the values is given in the header of the 
    file. At the end of the file there is also the number of objects tested, the
    number of objects detected, and the number of objects tracked.
  - The second (`SimulTrack.dat` in the example `Driver.in`) contains all the tracked
    objects. The meaning of the values is given in the header of the file.

    The parameters of the survey blocks determine which of 
    the detected objects are NOT tracked, based on the properties of the 
    non-tracked objects in the survey characterization used as input.

    NOTE: in both files, the color is <mag in survey filter> - <mag in
    reference "x" filter>, and the magnitude <mag> is given in the reference
    "x" filter.

---

### Model File Format  (aka lookup tables)
  The file describing the model must follow a strict format to be used with the
  routine GiMeObj provided in ReadModelFromFile.f. It MUST provide the
  following information to the program:

  - Epoch of elements
  - Colors of the objects
  - Orbital elements and absolute magnitude of the objects

  The absolut magnitude of an object is given in an arbitrary "x" band filter,
  which need not be in the list of known filters (g, r, i, z, u, B, V, R, I).
  Colors are given as "known_filter - x".

  The format of the model file is as follows (units are `au` and `degree`)


```
# Epoch of elements: JD = <Epoch_of_elements (JD)>
# Colors = <g-x> <r-x> <i-x> <z-x> <u-x> <V-x> <B-x> <R-x> <I-x> <?-x>
<a> <e> <inc> <Node> <Peri> <M> <H>
<a> <e> <inc> <Node> <Peri> <M> <H>
...
```

with

-    a     : semi-major axis [AU]
-    e     : eccentricity
-    inc   : inclination [degree]
-    Node  : longitude of node [degree]
-    Peri  : argument of peri-centre [degree]
-    M     : mean anomaly [degree]
-    H     : absolute magnitude [mag]

  The survey simulator expects orbital elements with respect to barycentric
  ecliptic reference frame, so the model must provides them in that reference
  frame.

  In this file, comments can be added by putting a "#" sign as the first
  character of a line. Any line starting with a "#" sign, except for those
  mentioned above, will be ignored.

---

## IMPLEMENTING YOUR OWN MODEL

  The workflow of the survey simulator driver is as follows:

    Loop (until told not to):
           call GiMeObj(arg_list_1)
           Check for GiMeObj failure; set exit if so
           call detection routine  Detos1(arg_list_2)
           Log detections
           Check for exit conditions
    Go back and loop

  `GiMeObj` is the critical routine that contains the model of the outer Solar
  System population(s) exposed to the Survey Simulator (and ultimately compared
  to real detections). On each call it must return one object, defined by the
  orbit, photometry, and lightcurve parameters described below (see also
  [API](#api)). It must accept a file name that tells the routine where to find
  needed parameters (if any) and a random number generator seed.

  When the driver calls `GiMeObj`, it expects one outer Solar System object
  defined by:

  1. An osculating **barycentric ecliptic J2000** orbit.
     The element set `(a, e, i, long_node, arg_peri, M, JD)` must be used, where
     the mean anomaly `M` is given at epoch `JD`. Propagation of the object
     position in `Detos1` is an unperturbed barycentric two-body problem.

  2. An absolute H magnitude and a set of colours in the major filters.
     H may be specified in any filter (here called `x`). Apparent magnitudes
     use the Bowell HG formalism, so Bowell's `G` is also supplied as `gb`.
     The colour array must at least cover every filter used by the survey
     blocks in the simulation. For example, if all blocks are in `g` and H is
     also in `g` (`x = g`), all colour terms may be zero. If imaging is in
     `g` or `r` only, one can supply `H_g` with `g-g = 0` (`color(1)`) and
     `r-g = -0.7` (`color(2)`). The band `x` need not be one of the nine
     predefined bands. `Detos1` returns magnitudes in the user-chosen `x` band.

  3. A lightcurve amplitude, period, and phase.
     Turn the lightcurve off by setting the amplitude to zero (the period must
     **not** be zero). Phase is the rotational phase at the orbital-element
     epoch (the JD at which `M` is given). Amplitude is peak-to-peak.

  The package provides two implementations of `GiMeObj`:

  - `InnerHotModel.f95` — analytical model that reads parameters from a file
    and generates objects as requested
  - `ReadModelFromFile.f95` — reads objects from a file, returns them one at a
    time, and signals when it reaches end of file

  To use your own model, put the routine in a file and create a symbolic link
  (or use `make Driver GIMEOBJ=MyGiMeObj`). `Driver.f95` includes `GiMeObj.f95`:

    `ln -s <whatever.f95> GiMeObj.f95`

  Fortran logical unit numbers **7 to 19** are reserved for `Driver` and
  `SurveySubs` and must not be used by `GiMeObj` or any routine you add to the
  driver. Use unit numbers from **20** upward.

  It is good practice that when first started, `GiMeObj` writes a file
  describing the model used, the version, and the date of the routine.

  Since this routine is called once for every object created, it should read
  required parameters on the first call and save them for later calls.

---

### API

The API (list of arguments, `arg_list_1` above) for `GiMeObj` is

    (filena, seed, a, e, inc, node, peri, M, epoch, h, color,
     gb, ph, period, amp, comment, nchar, ierr)

with:

#### INPUT
    filena: name of the file to be read in by GiMeObj the first time it
            is called (CH)
    seed  : Random number generator seed (I4)

#### OUTPUT
    a           : semimajor axis (R8)
    e           : eccentricity (R8)
    inc         : Inclination with respect to J2000 ecliptic [rad] (R8)
    node        : Longitude of ascending node [rad] (R8)
    peri        : Argument of perihelion [rad] (R8)
    M           : Mean anomaly [rad] (R8)
    epoch       : epoch for M (and rotational phase below), in Julian Day (R8)
    h           : absolute magnitude of object in band filter "x" (R8)
    color       : array of colors "y-x", where the index of "y" is as
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
    gb          : opposition surge factor G, Bowell formalism (R8)
    ph          : phase of lightcurve at epoch [rad] (R8)
    period      : period of lightcurve [day] (R8)   CANNOT SET TO ZERO
    amp         : peak-to-peak amplitude of lightcurve [mag] (R8)
                  CAN  SET TO ZERO
    comment     : user specified string containing whatever the user wants
                  (CHAR*100); can be empty.
    nchar       : number of characters in the comment string that should be
                  printed out in output files if the object is detected;
                  maximum of 100 (I4)
    ierr        : return code
                     0 : GiMeObj does not diagnose any errors, normal return
                         value
                   100 : end of model, exit after checking this object
                   -10 : could not get all orbital elements, skip object
                   -20 : something went grossly wrong, should quit

  Normally the Driver terminates the simulator (enough tracked detections), but
  these return codes also let `GiMeObj` tell the driver when to stop.

---

### Detos1

  Attempt to DETect 1 Outer Solar-system object. After `GiMeObj` returns an
  object, the driver asks `Detos1`: given the detection efficiencies and
  pointing history of all survey blocks, and allowing for probabilistic
  detection (especially for faint objects), is this object in the field
  coverage and detected by any block? If so, where and when, how bright, and
  was it tracked to a high-precision orbit?

  The list of arguments (`arg_list_2` above) for `Detos1` is

    (a, e, inc, node, peri, mt0, jday, hx, color, gb, ph, period, amp, surnam,
     seed, flag, ra, dec, d_ra, d_dec, r, delta, mi_int, m_rand, eff, isur, mt,
     jdayp, ic, surna, h_rand)

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
    m_int : Intrinsic apparent magnitude, in x-band (R8); from absolute
            magnitude, Sun/Earth distance, and phase angle (Bowell formalism);
            returned in the user-defined 'x' filter
    m_rand: Averaged randomized magnitude, in x-band (R8); from the intrinsic
            magnitude plus Gaussian noise from `<mag_error>` in the efficiency
            files; returned in the user-defined 'x' filter
    eff   : Efficiency of detection of object (function of mag and survey) (R8)
    isur  : Identification number of the survey block the object was in (I4)
    mt    : Mean anomaly at discovery [rad] (R8)
    jdayp : Time of discovery [JD] (R8)
    ic    : Index of color used for survey (I4)
    surna : Detection survey name, for information (CH10)
    h_rand: Absolute randomized magnitude (R8)

