U	= poly_dec poly_lib datadec ioutils effut elemutils rot DE_size DE_const DE_chrhdr DE_ephhdr DE_utils xvutils getsur numutils surveysub
OBJU	= $(addsuffix .o,${U})
SRCU	= $(addsuffix .f95,${U})
FPPU	= $(addsuffix .fpp,${U})
MODU	= $(addsuffix .mod,${U})

M   = modelutils GiMeObj
OBJM	= $(addsuffix .o,${M})
SRCM	= $(addsuffix .f95,${M})
FPPM	= $(addsuffix .fpp,${M})
MODM	= $(addsuffix .mod,${M})


A	= poly_dec datadec ioutils elemutils rot
OBJA	= $(addsuffix .o,${A})
SRCA	= $(addsuffix .f95,${A})
FPPA	= $(addsuffix .fpp,${A})
MODA	= $(addsuffix .mod,${A})


SSIM	= SurveySubsF95
GOBJ	= GiMeObjF95

FC	= gfortran
FPP	= gfortran -E -x f95-cpp-input -fPIC -DPYTHON

FFLAGS	= -O3 -fPIC -x f95-cpp-input
FFLAGP	= -O3
CFLAGS  = 
CPPFLAGS = 

TOPTARGETS := clean test
# use the GIMEOBJ from the command line or the default value of InnerHotModel
GIMEOBJ ?= InnerHotModel

$(TOPTARGETS): $(SUBDIRS)

.PHONY: $(SUBDIRS) 

$(SUBDIRS):
	$(MAKE) -C $@ $(MAKECMDGOALS)

all:	Driver $(SUBDIRS)

test: 
	$(MAKE) -s -C ../tests test LANGUAGE=F95

%.fpp :	%.f95
	$(FPP) $< -o $@

%.o : %.f95
	$(FC) $(FFLAGS) -c $<

.PHONY: clean

clean: rmtmp
	rm -f _$(SSIM)*.so _$(GOBJ)*.so  _$(IOBJ)*.so _$(ROBJ)*.so Driver GiMeObj.f95

rmtmp:
	rm -f *~ *.o *.mod *.fpp f90wrap*f90 SimulDetect.dat SimulTrack.dat
	rm -f $(IOBJ).py $(ROBJ).py $(SSIM).py $(GOBJ).py
	rm -f *.log LOG .f2py_f2cmap

.PHONY: link

link:
	\rm -f GiMeObj.f95
	\rm -f GiMeObj.o
	ln -s $(GIMEOBJ).f95 GiMeObj.f95

Driver: link Driver.f95 $(OBJU) ${SRCM} Makefile
	$(FC) $(FFLAGP) -o Driver $(OBJU) ${SRCM} Driver.f95

magerr: link magerr.f95 $(OBJU) ${SRCM} Makefile
	$(FC) $(FFLAGP) -o magerr $(OBJU) ${SRCM} magerr.f95
	rm GiMeObj.f95

Modules: _$(SSIM).so Makefile
	echo "Modules have been built"

_$(SSIM).so: $(FPPU) $(OBJU) Makefile
	\rm -f _$(SSIM).so
	f90wrap -k kind_map -m $(SSIM) $(FPPU) 1>f90wrap_SurveySubs.log 2>&1
	f2py-f90wrap --fcompiler=$(FC) -c -m _$(SSIM) $(OBJU) f90wrap_*.f90 1>f2py_SurveySubs.log 2>&1
	\rm f90wrap_*.f90
	\rm *.fpp *.mod *.o
	sed -e 's/import _$(SSIM)/from . import _$(SSIM)/' -i.bck $(SSIM).py
	rm $(SSIM).py.bck
	rm .f2py_f2cmap
	mv _$(SSIM)*.so _$(SSIM).so
