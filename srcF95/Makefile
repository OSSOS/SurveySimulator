OBJU	= poly_dec.o poly_lib.o datadec.o ioutils.o effut.o elemutils.o rot.o xvutils.o getsur.o numutils.o surveysub.o
SRCU	= poly_dec.f95 poly_lib.f95 datadec.f95 ioutils.f95 effut.f95 elemutils.f95 rot.f95 xvutils.f95 getsur.f95 numutils.f95 surveysub.f95

OBJM	= gimeobjut.o
SRCM	= gimeobjut.f95

FC	= gfortran

FFLAGS	= -O3

all:	test_read_pointings Driver

%.o : %.f95
	$(FC) $(FFLAGS) -c $<

clean:
	rm -f *~ *.o *.mod SimulDetect.dat SimulTrack.dat

test_read_pointings: test_read_pointings.f95 $(OBJU) Makefile
	$(FC) $(FFLAGS) -o test_read_pointings $(OBJU) test_read_pointings.f95

Driver: Driver.f95 $(OBJU) $(OBJM) Makefile
	$(FC) $(FFLAGS) -o Driver $(OBJU) $(OBJM) Driver.f95
