
# FORTRAN Usage Examples/Tests

This directory contains the input files to run the survey simulator and test the installation.
Results are compared with check files generated on our development machine as a unit test.

To do the test run the following commands:

`./run_tests.sh`

## The content of this directory: 

- README.md (this file)           # Describes directory content
- run_tests.sh                    # script to run the test suite
- test.sh                         # script to run the test suite
- Makefile                        # Execute the unit test `make test LANGUAGE=(F77|F95)` (run by run_tests.sh)
- Models                          # Input model files used for the tests
- Surveys                         # Survey desciption files used in the test, look here to build your own
- InnerHotModel.in                # Input file to run SurveySimulator from a set of parameters.
- InnerHotModel-check-F77.dat     # Result that should be produced by InnerBeltModel.in when running F77 version
- InnerHotModel-check-F95.dat     # Result that should be produced by InnerBeltModel.in when running F95 version
- ReadModelFromFile.in            # Input file to run SurveySimulator for an orbt database file.
- ReadModelFromFile-check-F77.dat # Result that should be produced by ReadModelFromFile.in when running F77 version
- ReadModelFromFile-check-F95.dat # Result that should be produced by ReadModelFromFile.in when running F95 version
