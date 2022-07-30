
--------------------------------------------------------------------------------

This directory contains the input files to run the survey simulator with a
parametric model (the inner classical Kuiper Belt here), and compare the
results to what was generated on our development machine, as a mean to ensure
that the package copied correctly and compilation in your environment went
well.

The content of this directory is:

lookup -------- README.example (this file)  Describes directory content
       \------- Driver.in      	     	    Input file for driver
        \------ InnerHot.in		    Parameters for the model in
         \                                  GiMeObj
          \---- SimulDetect-check.dat	    Result of reference simulation
           \--- test.sh			    Run the survey simulator and
                                            compare with reference simulation;
                                            Assumes the executable is in the
                                            parent directory,
                                            ../SurveySimulator
  
--------------------------------------------------------------------------------
