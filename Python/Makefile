# CFEPS Survey Simulator Makefile
#

### also works with ifort and f77
FC = /usr/bin/gfortran

Modules: SurveySubs.so GiMeObj.so Makefile
	echo "Modules have been built"

SurveySubs.so: SurveySubs.f Makefile
	f2py -c --f77exec=$(FC) --f77flags=-fPIC -m SurveySubs SurveySubs.f

GiMeObj.so: InnerBeltModel.f SurveySubs.so Makefile
	\rm -f GiMeObj.f
	ln -s InnerBeltModel.f GiMeObj.f
	f2py -c --f77exec=$(FC) --f77flags=-fPIC -m GiMeObj ./SurveySubs.so GiMeObj.f
	\rm GiMeObj.f

clean:
	\rm -f *.o core SurveySimulator

test: SurveySimulator
	./example/test.sh
