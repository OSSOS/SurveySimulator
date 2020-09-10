# CFEPS Survey Simulator Makefile
#
### also works with ifort and f77
FC = gfortran
SUBDIRS := fortran python
TOPTARGETS := all clean test

$(TOPTARGETS): $(SUBDIRS)

$(SUBDIRS):
	$(MAKE) -C $@ $(MAKECMDGOALS)

.PHONY: $(TOPTARGETS) $(SUBDIRS)


