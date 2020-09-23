# Makefile for the Ephemeris_FM F90 programs

TOP = ./
SRCDIR = $(TOP)
LIBDIR = $(TOP)
BINDIR = $(TOP)

FC=gfortran
FFLAGS=-u -O3 -funroll-loops -malign-double -march=native -ffree-form -ffree-line-length-none

.SUFFIXES:
.SUFFIXES: .o .f90

.PHONY: clean

SRC_F = \
	JMP_SSEphem_module.f90 \
	check_module_JMP.f90

OBJ_F = $(SRC_F:.f90=.o)
MOD_F = $(SRC_F:.f90=.mod)

%.o: %.f90
	$(FC) $(FFLAGS) -c $<

check_module_JMP: check_module_JMP.f90 JMP_SSEphem_module.o
	$(FC) $(FFLAGS) -o $@ $< JMP_SSEphem_module.o

clean:
	-rm *.o *.mod check_module_JMP
