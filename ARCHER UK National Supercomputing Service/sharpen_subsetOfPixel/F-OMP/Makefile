MF=	Makefile

#COMPILER
FC=	ftn
# On ARCHER, ftn is a wrapper for whichever Fortran compiler - GNU (gfortran), Intel (ifort), or Cray (crayftn) - has been chosen by loading the appropriate PrgEnv module prior to compilation.
CC=	cc
# Likewise, cc is a wrapper for whichever C compiler - GNU (gcc), Intel (icc), or Cray (craycc) has been loaded.

#For the Data Analytic Cluster use gfortran and gcc
#FC=	gfortran
#CC=	gcc

FFLAGS=	-g -DF_OPENMP_PRACTICAL
#For the Data Analytic Cluster, enable OpenMP
#FFLAGS=	-g -DF_OPENMP_PRACTICAL -fopenmp
CFLAGS= $(FFLAGS)
LFLAGS=

EXE=	sharpen

FSRC= \
	sharpen.f90	\
	dosharpen.f90 	\
	filter.f90	\
	fio.f90

CSRC=	utilities.c

#
# No need to edit below this line
#

.SUFFIXES:
.SUFFIXES: .f90 .c .o

OBJ=	$(FSRC:.f90=.o) $(CSRC:.c=.o)

.f90.o:
	$(FC) $(FFLAGS) -c $<

.c.o:
	$(CC) $(CFLAGS) -c $<

all:	$(EXE)

$(EXE):	$(OBJ)
	$(FC) $(FFLAGS) -o $@ $(OBJ) $(LFLAGS)

$(OBJ):	$(MF)

tar:
	tar cvf $(EXE).tar $(MF) $(SRC) fuzzy.pgm submit.sh

clean:
	rm -f $(OBJ) $(EXE) core
