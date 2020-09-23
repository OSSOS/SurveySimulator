
all: test

test: clean test_IBM test_RMFF

test_IBM: InnerBeltModel
	./test.sh InnerBeltModel

test_RMFF: ReadModelFromFile
	./test.sh ReadModelFromFile
	
InnerBeltModel: ../fortran/InnerBeltModel.f95
	$(MAKE) -C ../fortran Modules GIMEOBJ=InnerBeltModel

ReadModelFromFile: ../fortran/ReadModelFromFile.f95
	$(MAKE) -C ../fortran Modules GIMEOBJ=ReadModelFromFile

clean:
	\rm -f ModelUsed.dat SimulDetect.dat SimulTrack.dat LOG
