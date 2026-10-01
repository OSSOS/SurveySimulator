# Things that are needed fort the CLI (User interface)
U	= modelutils GiMeObj 
OBJU	= $(addsuffix .o,${U})
SRCU	= $(addsuffix .f95,${U})
FPPU	= $(addsuffix .fpp,${U})
MODU	= $(addsuffix .mod,${U})

# Things that should end up with python module access
M   = datadec ioutils surveysub
OBJM	= $(addsuffix .o,${M})
SRCM	= $(addsuffix .f95,${M})
FPPM	= $(addsuffix .fpp,${M})
MODM	= $(addsuffix .mod,${M})
F90M    = $(addsuffix .f90,${M})
F90_WRAP_M    = $(addprefix f90wrap_,${F90M})

# Things that should be in the library
A   = parameters poly_dec datadec common_data debug poly_lib ioutils elemutils rot effut xvutils numutils getsur surveysub
OBJA	= $(addsuffix .o,${A})
SRCA	= $(addsuffix .f95,${A})
FPPA	= $(addsuffix .fpp,${A})
MODA	= $(addsuffix .mod,${A})

MODULE	?= ossssimlib
GOBJ	= GiMeObjF95


# FC is read by NumPy's meson f2py backend (Python>=3.12); --fcompiler is ignored there.
# Make's built-in default is FC=f77; ?= would not override it.
ifeq ($(origin FC),default)
FC := gfortran
endif
FC	?= gfortran
FPP	= $(FC) -E -x f95-cpp-input -fPIC -DPYTHON
# FPP_F90FLAGS = -x f95-cpp-input -fPIC


FFLAGS	= -O3 -fPIC -x f95-cpp-input
FFLAGP	= -O3
CFLAGS  = 
CPPFLAGS = 
LIBTOOL = ar src
# Temporary tree for numpy.f2py's meson backend (not a project-level Meson build).
BUILD_DIR = f2py_build
KIND_MAP = kind_map
F2PYFLAGS = --build-dir $(BUILD_DIR)

# NumPy f2py meson.build lists -lsrc before -L…; GNU ld then fails to find libsrc.
# LIBRARY_PATH resolves the archive regardless of flag order. LDFLAGS/FCFLAGS help
# meson and the Fortran compiler find libsrc.a and .mod files from bbdir.
export LIBRARY_PATH := $(CURDIR)$(if $(LIBRARY_PATH),:$(LIBRARY_PATH),)
export LDFLAGS := -L$(CURDIR)$(if $(LDFLAGS), $(LDFLAGS),)
export FCFLAGS := -I$(CURDIR)$(if $(FCFLAGS), $(FCFLAGS),)


.PHONY: all clean

# TOPTARGETS := clean
# use the GIMEOBJ from the command line or the default value of InnerHotModel
GIMEOBJ ?= InnerHotModel

# $(TOPTARGETS): $(SSIM)

# .PHONY: $(SUBDIRS) 

# $(SUBDIRS):
# 	$(MAKE) -C $@ $(MAKECMDGOALS)

all: $(MODULE)

test: 
	$(MAKE) -s -C tests test LANGUAGE=F95

%.fpp :	%.f95
	$(FPP) $< -o $@

%.o : %.f95
	$(FC) $(FFLAGS) -c $<

.PHONY: clean

clean: rmtmp
	rm -f _$(MODULE)*.so _$(GOBJ)*.so  _$(IOBJ)*.so _$(ROBJ)*.so Driver GiMeObj.f95
	rm -rf $(BUILD_DIR) libsrc.a

rmtmp:
	rm -f *~ *.o *.mod *.fpp f90wrap*f90 SimulDetect.dat SimulTrack.dat
	rm -f $(IOBJ).py $(ROBJ).py $(MODULE).py $(GOBJ).py _$(MODULE)*.so
	rm -f *.log LOG .f2py_f2cmap

.PHONY: link

link:
	\rm -f GiMeObj.f95
	\rm -f GiMeObj.o
	ln -s $(GIMEOBJ).f95 GiMeObj.f95

Driver: link Driver.f95 $(OBJA) $(OBJU) Makefile GiMeObj.f95
	$(FC) $(FFLAGP) -o Driver $(OBJA) $(OBJU) Driver.f95

magerr: link magerr.f95 $(OBJU) ${SRCM} Makefile
	$(FC) $(FFLAGP) -o magerr $(OBJU) ${SRCM} magerr.f95
	rm GiMeObj.f95

libsrc.a: ${OBJA}
	${LIBTOOL} $@ $?

Modules: _$(MODULE).so Makefile 
	echo "Modules have been built"

$(MODULE): Makefile $(FPPM) libsrc.a
	# Force fresh .mod/.o visibility for f90wrap; stale debug.mod caused CANFAR clashes.
	# Wipe any prior f2py/meson tree so we never reuse deleted /tmp/pip-build-env-* paths.
	rm -rf $(BUILD_DIR)
	f90wrap -m $(MODULE) $(FPPM) -k $(KIND_MAP) -v -P > f90wrap.log 2>&1
	mkdir -p $(BUILD_DIR)
	# .mod files must be visible from f2py's bbdir (sibling of these symlinks).
	ln -s -f $(CURDIR)/*.mod $(BUILD_DIR)/
	# Surface the real compiler/linker error; pip only shows make's exit status otherwise.
	# Prefer f2py_f90wrap.log (ninja) over meson-log.txt (configure-only, often stale).
	f2py-f90wrap $(F2PYFLAGS) -c -m _$(MODULE) $(F90_WRAP_M) -L$(CURDIR) -lsrc > f2py_f90wrap.log 2>&1 \
		|| { \
			echo "==== f2py_f90wrap.log (last 80 lines) ===="; \
			tail -n 80 f2py_f90wrap.log; \
			if [ -d $(BUILD_DIR)/bbdir ]; then \
				echo "==== meson compile -C $(BUILD_DIR)/bbdir -v ===="; \
				meson compile -C $(BUILD_DIR)/bbdir -v || true; \
			fi; \
			if [ -f $(BUILD_DIR)/bbdir/meson-logs/meson-log.txt ]; then \
				echo "==== meson-log.txt (last 40 lines) ===="; \
				tail -n 40 $(BUILD_DIR)/bbdir/meson-logs/meson-log.txt; \
			fi; \
			exit 1; \
		}
	\rm f90wrap*.f90 *.fpp *.o *.mod .f2py_f2cmap libsrc.a
	\rm -r $(BUILD_DIR) 

