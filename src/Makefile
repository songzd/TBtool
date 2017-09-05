.SUFFIXES: .x .o .f90 .F

FC = ifort 
MPI = ifort
FLAGS = -mkl -g -warn -check bounds 
#LIBS = -L/home/soft/intel/mkl/lib/intel64 \
       -lmkl_intel_lp64 -lmkl_sequential -lmkl_core

SUFFIX=.f90
CPP=fpp -f_com=no -free -w0 $*.F $*$(SUFFIX)


MODULE = typcnst.o my_math.o pntgrp.o reprsnt.o band.o tghtbd.o symm.o slicewann.o
EXE =  TBtool

default : $(EXE)

$(EXE):%: $(MODULE) %.o
	$(FC) $(FLAGS) $(LIB) -o  $@ $^

%.o:%.F
	$(CPP)
	$(FC) -c $(FLAGS) $*.f90


clean:
	rm -f *.o *.mod *.f90
