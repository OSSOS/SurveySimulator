U	= poly_dec poly_lib datadec ioutils effut elemutils rot xvutils getsur numutils surveysub
OBJU	= $(addsuffix .o,${U})
SRCU	= $(addsuffix .f95,${U})
FPPU	= $(addsuffix .fpp,${U})
MODU	= $(addsuffix .mod,${U})

M	= modelutils gimeobjut
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
FPP	= gfortran -E -x f95-cpp-input -fPIC

FFLAGS	= -O3 -fPIC -x f95-cpp-input
FFLAGP	= -O3

all:	test_read_pointings Driver Modules

%.fpp :	%.f95
	$(FPP) $< -o $@

%.o : %.f95
	$(FC) $(FFLAGS) -c $<

clean:
	rm -f *~ *.o *.mod *.fpp f90wrap*f90 $(SSIM).py $(GOBJ).py _$(SSIM)*.so _$(GOBJ)*.so SimulDetect.dat SimulTrack.dat

test_read_pointings: test_read_pointings.f95 $(OBJU) Makefile
	$(FC) $(FFLAGP) -o test_read_pointings $(OBJU) test_read_pointings.f95

Driver: Driver.f95 $(SRCM) $(OBJU) $(OBJM) Makefile
	$(FC) $(FFLAGP) -o Driver $(OBJU) $(OBJM) Driver.f95

Modules: _$(SSIM).so _$(GOBJ).so Makefile
	echo "Modules have been built"

_$(SSIM).so: $(FPPU) $(OBJU) Makefile
	\rm -f _$(SSIM).so
	f90wrap -k kind_map -m $(SSIM) $(FPPU) 1>f90wrap_SurveySubs.log 2>&1
	f2py-f90wrap --fcompiler=$(FC) -c -m _$(SSIM) $(OBJU) f90wrap_*.f90 1>f2py_SurveySubs.log 2>&1
	\rm f90wrap_*.f90
	mv _$(SSIM)*.so _$(SSIM).so

_$(GOBJ).so: $(FPPM) $(OBJM) $(FPPA) $(OBJA) _$(SSIM).so Makefile
	\rm -f _$(GOBJ).so
	f90wrap -k kind_map -m $(GOBJ) $(FPPM) $(FPPA) 1>f90wrap_GiMeObj.log 2>&1
	f2py-f90wrap --fcompiler=$(FC) -c -m _$(GOBJ) $(OBJA) $(OBJM) f90wrap_*.f90 1>f2py_GiMeObj.log 2>&1
	\rm f90wrap_*.f90
	mv _$(GOBJ)*.so _$(GOBJ).so
