# CFEPS Survey Simulator Makefile
#
### also works with ifort and f77
FC = gfortran

SUBDIRS := example
TOPTARGETS := all clean test

$(TOPTARGETS): $(SUBDIRS)

$(SUBDIRS):
	$(MAKE) -C $@ $(MAKECMDGOALS)

.PHONY: $(TOPTARGETS) $(SUBDIRS)

%.o : %.f
	$(FC) $(FFLAGS) -c $<

all: Driver

link:
	\rm -f GiMeObj.f
	ln -s $(GIMEOBJ).f GiMeObj.f

Driver: link Driver.f SurveySubs.f GiMeObj.o Makefile
	$(FC) -O3 GiMeObj.f SurveySubs.f Driver.f -o Driver

clean:
	\rm -f *.o *.so core Driver *.log LOG GiMeObj.f

# Python Module parts

Modules: SurveySubs.so ModelUtils.so GiMeObj.so Makefile
	echo "Modules have been built"

SurveySubs.so: SurveySubs.f Makefile
	f2py -c --f77exec=$(FC) --f77flags=-fPIC --f90flags=-fPIC -m SurveySubs SurveySubs.f 1>f2py_SurveySubs.log 2>&1
	mv SurveySubs.*.so SurveySubs.so

ModelUtils.so: ModelUtils.f Makefile
	f2py -c --f77exec=$(FC) --f77flags=-fPIC --f90flags=-fPIC -m ModelUtils ModelUtils.f 1>f2py_ModelUtils.log 2>&1
	mv ModelUtils.*.so ModelUtils.so

GiMeObj.so: link GiMeObj.f Makefile
	f2py -c --f77exec=$(FC) --f77flags=-fPIC -m GiMeObj GiMeObj.f 1>f2py_GiMeObj.log 2>&1
	\rm GiMeObj.f
	mv GiMeObj.*.so GiMeObj.so

