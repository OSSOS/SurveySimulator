
all: test

test: clean test_IBM test_RMFF

test_IBM: InnerBeltModel
	./test.sh InnerBeltModel

test_RMFF: ReadModelFromFile
	./test.sh ReadModelFromFile

F2PY = f90wrap

K := $(if $(shell which $(F2PY)),0,$(error "No $(F2PY) in PATH. Please install f90wrap"))
	
InnerBeltModel: ../fortran/InnerBeltModel.f95 
	$(MAKE) -C ../fortran Modules GIMEOBJ=InnerBeltModel

ReadModelFromFile: ../fortran/ReadModelFromFile.f95
	$(MAKE) -C ../fortran Modules GIMEOBJ=ReadModelFromFile

clean:
	\rm -f ModelUsed.dat SimulDetect.dat SimulTrack.dat LOG test .f2py_f2cmap
	\rm -rf __pycache__
