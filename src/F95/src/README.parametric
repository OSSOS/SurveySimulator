
--------------------------------------------------------------------------------

This directory contains the source code for a GiMeObj routine that generates
objects according to some parametric presciption of the outer Solar System
population being modeled, along with a Makefile to generate the executable 
and an example subdirectory.

The contents of this directory is:

parametric ----- README.parametric (this file)  Describes "GiMeObj" routine
 \         \---- InnerHotModel.f                Example "GiMeObj" routine to
 |          \                                   generate objects via parametric
 |           \                                  prescription of inner belt
 |            \- Makefile                       Simplistic Makefile to build
 |                                              executable
 |
 \----example                                   Contains files to run an
                                                example and check validity of
                                                results

--------------------------------------------------------------------------------

COMPILING
  The simplest way to create the executable for the survey simulator is 
  to type the following command at the prompt:

    make

  This will create an executable file called "SurveySimulator".  If you do
  not have the gfortran compiler installed you will need to replace the
  compiler command in the Makefile with the fortran compiler you use.

  One can compile the program directly. Because Driver.f 'includes' a
  file 'GiMeObj.f' containing the model definition, one must first create 
  a symbolic link pointing to the actual file "InnerHotModel.f":

    ln -s InnerHotModel.f GiMeObj.f
    gfortran -O3 -o SurveySimulator Driver.f

  Note: Driver.f SurveySubs.f in this directory are symbolically linked 
  to the actual source codes in the src subdirectory.

RUNNING
  The driver program provided with this package reads six parameters from the
  standard input:

  - the seed for the random number generator (integer)
  - a number <n> to control how long we run:
      <n> > 0: maximum number of simulated tracked detections
      <n> < 0: -maximum number of iterations (i.e. the number of calls to
               GiMeObj)
      <n> = 0: run until the 'model' decides to stop
  - the name of the directory containing the characterized survey blocks
  - the name of the model input file (a file GiMeObj will read in)
  - the name of the output file where the detected objects will be listed
  - the name of the output file where the detected and tracked objects will be
    listed


  An easy way of feeding these values into the driver is to write
  them in a file, say "Driver.in", one line for each value, and execute the
  program with redirection of standard input:

    SurveySimulator < Driver.in

  An example of this Driver.in file is in the parametric/example subdirectory.

  Execution generates two output files with the user specified names:
  - The first (called SimulDetect.dat in the example cases given) contains all
    detected objects. The meaning of the values is given in the header of the 
    file. At the end of the file there is also the number of objects tested, the
    number of objects detected, and the number of objects tracked.
  - The second (called SimulTrack.dat in the examples) contains all the tracked
    objects. The meaning of the values is given in the header of the file.

    The parameters of the survey blocks determine which (small) fraction of
    the detected objects are NOT tracked, based on the properties of the 
    non-tracked objects in the true survey.

    NOTE: in both files, the color is <mag in survey filter> - <mag in
    reference "x" filter>, and the magnitude <mag> is given in the reference
    "x" filter.

--------------------------------------------------------------------------------

IMPLEMENTING YOUR OWN MODEL

  As described in ../src/README.src, the workflow of the survey simulator
  driver is as follows:

    Loop (until told not to):
           call GiMeObj(arg_list_1)
           Check for GiMeObj failure; set exit if so
           call detection routine  Detos1(arg_list_2)
           Log detections
           Check for exit conditions
    Go back and loop

  The GiMeObj routine is in charge of providing a single new object at each
  call, an object being defined by (see below) its orbital elements, the 
  absolute magnitude of the object in some band filter "x", the colors of the 
  object, the opposition surge effect parameter (G in Bowell's formalism) and 
  lightcurve parameters (period, peak-to-peak amplitude and phase at epoch). 
  GiMeObj must accept a file name as input that tells the routine where to 
  find the needed parameters (if any) and also a random number generator seed. 
  The package provides two different implementations of the GiMeObj. 
   -The one in InnerHotModel.f shows an example of an analytical model that 
    reads its parameters from a file and then generates objects as requested.
   -The one in FromFileModel.f will read objects from a file, return them 
    one at a time and signal when it has reach the end of the file. 

  One can easily create one's own GiMeObj routine to replace the ones provided
  with the package. The Driver.f program uses an 'include' for the file
  'GiMeObj.f' containing the model definition. The suggested way to use
  this feature with one's own code is to have one's GiMeObj routine in a file
  <whatever.f> and create a symobolic link:

    ln -s <whatever.f> GiMeObj.f

  The model subroutines GiMeObj can access files using Fortran logical unit
  numbers from 20 upward. This range in reseved for them and won't be used by
  the drivers nor SurveySubs routines.

  It is good practice that when first started, the GiMeObj routine writes a
  file describing the model used, the version and the date of the routine.

  Since this routine is called once for every object created, it needs to get
  all the required parameters once when it is called the first time, then save
  these values for future use.

  The survey simulator expects orbital elements with respect to barycentric
  ecliptic reference frame, so the model must provides them in that reference
  frame.

--------------------------------------------------------------------------------

The API (list of arguments, arg_list_1 above) for GiMeObj is

    (filena, seed, a, e, inc, node, peri, M, epoch, h, color,
     gb, ph, period, amp, comment, nchar, ierr)

with:

INPUT
    filena: name of the file to be read in by GiMeObj the first time it
            is called (CH)
    seed  : Random number generator seed (I4)

OUTPUT
    a           : semimajor axis (R8)
    e           : eccentricity (R8)
    inc         : Inclination with respect to J2000 ecliptic [rad] (R8)
    node        : Longitude of ascending node [rad] (R8)
    peri        : Argument of perihelion [rad] (R8)
    M           : Mean anomaly [rad] (R8)
    epoch       : epoch for M (and rotational phase below), in Julian Day (R8)
    h           : absolute magnitude of object in band filter "x" (R8)
    color       :varray of colors "y-x", where the index of "y" is as
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
    comment  	: user specified string containing whatever the user wants
		  (CHAR*100); can be empty.
    nchar	: number of characters in the comment string that should be
    		  printed out in output files if the object is detected;
                  maximum of 100 (I4)
    ierr        : return code
                     0 : GiMeObj does not diagnose any errors, normal return
                         value
                   100 : end of model, exit after checking this object
                   -10 : could not get all orbital elements, skip object
                   -20 : something went grossly wrong, should quit


--------------------------------------------------------------------------------
