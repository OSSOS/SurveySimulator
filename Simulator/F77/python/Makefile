
all: test

test: clean test_IBM 

test_IBM: InnerHotModel
	./test.sh InnerModelModel

Res32Model: ../fortran/Res32Model.f
	$(MAKE) -C ../fortran Modules GIMEOBJ=Res32Model

clean:
	\rm -f ModelUsed.dat SimulDetect.dat SimulTrack.dat LOG ModelUsed
