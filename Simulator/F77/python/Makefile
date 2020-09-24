
%: %.sh

all: test

test: clean test_IBM  test_RMFF test_R32M

test_IBM: InnerHotModel
	./test.sh InnerHotModel

test_RMFF: ReadModelFromFile
	./test.sh ReadModelFromFile

test_R32M: Res32Model
	./test.sh Res32Model

InnerHotModel: ../fortran/InnerHotModel.f
	$(MAKE) -C ../fortran Modules GIMEOBJ=InnerHotModel

ReadModelFromFile: ../fortran/ReadModelFromFile.f
	$(MAKE) -C ../fortran Modules GIMEOBJ=ReadModelFromFile

Res32Model: ../fortran/Res32Model.f
	$(MAKE) -C ../fortran Modules GIMEOBJ=Res32Model

clean:
	\rm -f ModelUsed.dat SimulDetect.dat SimulTrack.dat LOG ModelUsed
