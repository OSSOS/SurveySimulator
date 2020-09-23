
--------------------------------------------------------------------------------

This directory contains the source code for GiMeObj routine that reads a model
from a file (lookup table), along with a Makefile to generate the executable
and an example subdirectory.

The content of this directory is:

lookup -------- README.lookup (this file)  Describes format of model files
 \     \------- ReadModelFromFile.f        Defines "GiMeObj" routine to read a
 |      \                                  model from a file (lookup table)
 |       \----- Makefile                   Simplistic Makefile to build
 |                                         executable 
 |
 \----example                              Contains files to run an example and
                                           check validity of results

--------------------------------------------------------------------------------

COMPILING
  The simplest way to create the executable for the survey simulator is to type
  the following command at the prompt:

    make

  This will create an executable file called "SurveySimulator".

  One can also easily compile the program directly. Since the driver includes a
  file 'GiMeObj.f' containing the definition of the model, one must first
  create a symbolic link pointing to the actual file "ReadModelFromFile.f":

    ln -s ReadModelFromFile.f GiMeObj.f
    gfortran -O3 -o SurveySimulator Driver.f

RUNNING
  The driver program provided with this package reads six parameters from the
  standard input:

  - the seed for the random number generator (integer)
  - a number <n> to control how long we run:
      <n> > 0: maximum number of simulated tracked detections
      <n> < 0: -maximum number of iterations (i.e. the number of calls to
               GiMeObj)
      <n> = 0: run until the 'model' decides to stop
  - the name of the directory containing the survey characterization.
  - the name of the model input file (a file GiMeObj will read in)
  - the name of the output file where the detected objects will be listed
  - the name of the output file where the detected and tracked objects will be
    listed

  An easy and natural way of feeding these values into the driver is to write
  them in a file, say "Driver.in", one line for each value, and execute the
  program with redirection of standard input:

    SurveySimulator < Driver.in

  This will generate two output files with names specified by the user:
  - The first (called SimulDetect.dat in the example cases given) contains all  
    detected objects. The meaning of the values is given in the header of the 
    file. At the end of the file there is also the number of objects tested, 
    the number of objects detected, and the number of objects tracked.
  - The second (called SimulTrack.dat in the examples) contains all the tracked 
    objects. The meaning of the values is given in the header of the file.

    The parameters of the survey blocks determine which (small) fraction of
    the detected objects are NOT tracked, based on the properties of the 
    non-tracked objects in the true survey.

    NOTE: in both files, the color is <mag in survey filter> - <mag in
    reference "x" filter>, and the magnitude <mag> is given in the reference
    "x" filter.

--------------------------------------------------------------------------------

FORMAT OF MODEL FILES (aka lookup tables)
  The file describing the model must follow a strict format to be used with the
  routine GiMeObj provided in ReadModelFromFile.f. It MUST provide the
  following informations to the program:

  - Epoch of elements
  - Colors of the objects
  - Orbital elements and absolute magnitude of the objects

  The absolut magnitude of an object is given in an arbitrary "x" band filter,
  which need not be in the list of known filters (g, r, i, z, u, B, V, R, I).
  Colors are given as "known_filter - x".

  The format of the model file is as follows:

# Epoch of elements: JD = <Epoch_of_elements (JD)>
# Colors = <g-x> <r-x> <i-x> <z-x> <u-x> <V-x> <B-x> <R-x> <I-x> <?-x>
<a> <e> <inc> <Node> <Peri> <M> <H>
<a> <e> <inc> <Node> <Peri> <M> <H>
...

with
    a     : semi-major axis [AU]
    e     : eccentricity
    inc   : inlination [degree]
    Node  : longitude of node [degree]
    Peri  : argument of pericenter [degree]
    M     : mean anomaly [degree]
    H     : absolute magnitude [mag]

  The survey simulator expects orbital elements with respect to barycentric
  ecliptic reference frame, so the model must provides them in that reference
  frame.

  In this file, comments can be added by putting a "#" sign as the first
  character of a line. Any line starting with a "#" sign, except for those
  mentionned above, will be ignored.

--------------------------------------------------------------------------------
