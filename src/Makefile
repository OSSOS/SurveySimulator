# CFEPS Survey Simulator Makefile
#
### also works with ifort and f77
FC = gfortran

TOPTARGETS := all clean

SUBDIRS := test

$(TOPTARGETS): $(SUBDIRS)

$(SUBDIRS):
	$(MAKE) -C $@ $(MAKECMDGOALS)

.PHONY: $(TOPTARGETS) $(SUBDIRS)

all: Parametric ReadModelFromFile

Parametric: Driver.f SurveySubs.f InnerHotModel.f Makefile
	$(FC) -O3 InnerHotModel.f SurveySubs.f Driver.f -o Parametric

ReadModelFromFile: Driver.f SurveySubs.f ReadModelFromFile.f Makefile
	$(FC) -O3 ReadModelFromFile.f SurveySubs.f Driver.f -o ReadModelFromFile

clean:
	\rm -f *.o core Parametric ReadModelFromFile

