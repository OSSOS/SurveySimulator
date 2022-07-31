
# FORTRAN 95 Usage Examples and Tests

This directory contains the input files to run the survey simulator using input models. 
Results can be compared with check files generated on our development machine as a unit test.

To do the test run the following commands:

`make test`

## The content of this directory: 

README.md (this file)        # Describes directory content
test.sh                      # script to run the test suite
Makefile                     # Execute the unit test `make test` (run by test.sh)
InnerHotModel.in            # Input file to run SurveySimulator from a set of parameters.
InnerHotModel-check.dat     # Result that should be produced by InnerBeltModel.in 
ReadModelFromFile.in         # Input file to run SurveySimulator for an orbt database file.
ReadModelFromFile-check.dat  # Result that should be produced by ReadModelFromFile.in
