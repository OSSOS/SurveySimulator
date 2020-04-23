# CFEPS Survey Simulator Makefile
#

### also works with ifort and f77
FC = /usr/bin/gfortran

Modules: SurveySubs.so GiMeObj.so ModelUtils.so Makefile
	echo "Modules have been built"

SurveySubs.so: SurveySubs.f Makefile
	f2py -c --f77exec=$(FC) --f77flags=-fPIC --f90flags=-fPIC -m SurveySubs SurveySubs.f 1>f2py_SurveySubs.log 2>&1
	mv SurveySubs.*.so SurveySubs.so

ModelUtils.so: ModelUtils.f Makefile
	f2py -c --f77exec=$(FC) --f77flags=-fPIC --f90flags=-fPIC -m ModelUtils ./SurveySubs.so ModelUtils.f 1>f2py_ModelUtils.log 2>&1
	mv ModelUtils.*.so ModelUtils.so

GiMeObj.so: Res32Model.f SurveySubs.so Makefile
	\rm -f GiMeObj.f
	ln -s Res32Model.f GiMeObj.f
	f2py -c --f77exec=$(FC) --f77flags=-fPIC -m GiMeObj ./SurveySubs.so GiMeObj.f 1>f2py_GiMeObj.log 2>&1
	\rm GiMeObj.f
	mv GiMeObj.*.so GiMeObj.so

clean:
	\rm -f *.o core SurveySimulator

test: Modules
	./test.sh
