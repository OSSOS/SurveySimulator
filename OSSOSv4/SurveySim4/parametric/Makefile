# CFEPS Survey Simulator Makefile
#

### also works with ifort and f77
FC = gfortran

SurveySimulator: ../src/Driver.f ../src/SurveySubs.f InnerHotModel.f Makefile
	\rm -f GiMeObj.f
	ln -s InnerHotModel.f GiMeObj.f
	$(FC) -O3 -o SurveySimulator Driver.f
	\rm GiMeObj.f

clean:
	\rm -f *.o core SurveySimulator

test: SurveySimulator
	./example/test.sh
