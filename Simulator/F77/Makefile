# CFEPS Survey Simulator Makefile
#
### also works with ifort and f77
FC = gfortran
SUBDIRS := src
TOPTARGETS := all clean

$(TOPTARGETS): $(SUBDIRS)

$(SUBDIRS):
	$(MAKE) -C $@ $(MAKECMDGOALS)

.PHONY: $(TOPTARGETS) $(SUBDIRS)

all: $(SUBDIRS)

Modules: $(SUBDIRS)

test: $(SUBDIRS)
