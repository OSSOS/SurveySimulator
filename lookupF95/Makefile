# OSSOS Survey Simulator Makefile
#

all:	SurveySimulator

clean:
	\rm -f *.o *.mod core SurveySimulator

SurveySimulator:
	./CompileDriver.sh

test: SurveySimulator
	./example/test.sh
