# CFEPS Survey Simulator Makefile
#
### also works with ifort and f77
FC = gfortran

SUBDIRS := example
TOPTARGETS := all clean test
GIMEOBJ ?= InnerHotModel

#$(TOPTARGETS): $(SUBDIRS)

$(SUBDIRS):
	$(MAKE) -C $@ $(MAKECMDGOALS)

.PHONY: $(TOPTARGETS) $(SUBDIRS)

all: Driver

%.o : %.f
	$(FC) $(FFLAGS) -c $<

.PHONY: link

link:
	\rm -f GiMeObj.f
	\rm -f GiMeObj.o
	ln -s $(GIMEOBJ).f GiMeObj.f

Driver: link Driver.f SurveySubs.f GiMeObj.f Makefile
	$(FC) -O3 GiMeObj.f SurveySubs.f Driver.f -o Driver

.PHONY: clean

clean:
	\rm -f *.o *.so core Driver *.log LOG GiMeObj.f

# Python Module parts

Modules: SurveySubs.so ModelUtils.so GiMeObj.so Makefile
	echo "Modules have been built"

SurveySubs.so: SurveySubs.f Makefile
	f2py -c --f77exec=$(FC) --f77flags=-fPIC --f90flags=-fPIC -m SurveySubs SurveySubs.f 1>f2py_SurveySubs.log 2>&1
	mv SurveySubs.*.so SurveySubs.so

ModelUtils.so: ModelUtils.f SurveySubs.so Makefile
	f2py -c --f77exec=$(FC) --f77flags=-fPIC --f90flags=-fPIC -m ModelUtils SurveySubs.so ModelUtils.f 1>f2py_ModelUtils.log 2>&1
	mv ModelUtils.*.so ModelUtils.so

GiMeObj.so: link $(GIMEOBJ).f SurveySubs.so Makefile
	f2py -c --f77exec=$(FC) --f77flags=-fPIC -m GiMeObj SurveySubs.so GiMeObj.f 1>f2py_GiMeObj.log 2>&1
	\rm GiMeObj.f
	mv GiMeObj.*.so GiMeObj.so

