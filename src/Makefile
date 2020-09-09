# CFEPS Survey Simulator Makefile
#
### also works with ifort and f77
FC = gfortran
BASEDIR = $(CURDIR)

TOPTARGETS := all clean

SUBDIRS := F95 F77

$(TOPTARGETS): $(SUBDIRS)

$(SUBDIRS):
	$(MAKE) -C $@ $(MAKECMDGOALS)

.PHONY: $(TOPTARGETS) $(SUBDIRS)

all: $(SUBDIRS)

clean: $(SUBDIRS)

