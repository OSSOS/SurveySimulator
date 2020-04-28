
--------------------------------------------------------------------------------

This directory contains the common source codes for the Survey Simulator
(Driver.py and SurveySubs.f) and an example model (InnerBeltModel.f).
 

The content of this directory is:

src 
  \---------- README.python (this file)  Describes the architecture of the
   \                                     driver and GiMeOjb routines
    \-------- Driver.python              The driver program that rules them all
     \------- SurveySubs.f               A set of subroutines and functions that
      \                                  determine if an object was found, and
       \                                 also provide interesting utilities
        \---- InnerBeltModel.f           Example "GiMeObj" routine to generate
         \                               objects via parametric prescription
          \                              of inner belt
           \- Makefile                   Simplistic Makefile to build executable

--------------------------------------------------------------------------------

The workflow of the survey simulator driver is as follows:

    Loop (until told not to):
           call GiMeObj(arg_list_1)
           Check for GiMeObj failure; set exit if so
           call detection routine  Detos1(arg_list_2)
           Log detections
           Check for exit conditions
    Go back and loop

The GiMeObj routine (whose content should be determined by the user) returns 
a single new object at each call.

The Detos1 routine determines if the proposed object would have been detected
and possibly tracked (meaning the orbit of a detected object is determined to 
high precision) by the survey(s) the user has chosen to simulate. It then
decides if the object is characterized, i.e. is brighter than the
characterization limit for the block and rate of motion it was discovered
in. It is the characterized tracked simulated detection list which must be
compared to the true detections.

As currently written, Driver.py imports the module GiMeObj.so which contains
the definition of the model. The preferred way to use this is to have one's
GiMeObj routine in a file <whatever.f> and create a symobolic link to
GiMeObj.f:

    ln -s <whatever.f> GiMeObj.f

With this link, you can compile the model and SurveySubs.f and generate the
modules that are used by Driver.py.


COMPILING
  The simplest way to create the executable for the survey simulator is 
  to type the following command at the prompt:

    make

  This will create the two modules GiMeObj.so and SurveySubs.so imported by
  Driver.py. If you do not have the gfortran compiler installed you will need
  to replace the compiler command in Makefile with the fortran compiler you
  use.

  One can create the modules directly:

    ln -s InnerHotModel.f GiMeObj.f
    f2py -c --f77exec=/usr/bin/gfortran --f77flags=-fPIC -m SurveySubs SurveySubs.f
    f2py -c --f77exec=/usr/bin/gfortran --f77flags=-fPIC -m GiMeObj ./SurveySubs.so GiMeObj.f


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

  An easy way of feeding these values into the driver is to write them in a
  file, say "Driver.in", one line for each value, and execute the program with
  redirection of standard input:

    Driver.py < Driver.in

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
    reference "x" filter>, the intrinsic magnitudes <m_int> and <H_int> are
    given in the reference "x" filter, and the apparent magnitudes <m_rand> and
    <H_rand> are given in the detection filter.

--------------------------------------------------------------------------------
