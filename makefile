COMPILER=gfortran
#COMPILER = gfortran 		#compile with gfortran
#COMPILER = ifort 		#compile with ifort

FLAGS=-O4
	# optimization flags
	#FLAGS = -fbounds-check  #flag to check if array boundaries are respected\
	
ifeq ($(COMPILER),gfortran)	# gfortran compiler, standard LAPACK and BLAS linking
        LIB = -L/usr/lib/ -llapack -lblas
        # provide the appropriate path for LAPACK and BLAS libraries, they differ between distributions
        
else ifeq ($(COMPILER),ifort)	# ifort compiler, linking to mkl libraries
        MKLPATH = /share/apps/intel/mkl/lib/intel64
        MKLINCLUDE = /share/apps/intel/mkl/include
        LIB = -L$(MKLPATH) -I$(MKLINCLUDE) -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lfftw3xf_intel
        # provide the appropriate paths and flags to the MKL library, they differe for different versions of IFORT
        # the above choice works for IFORT version 14.0.2 20140120 
endif

main: mps_grand.x

mps_grand.x: open_dmrg.o bose_hubbard.o eigensolver.o eigensolver_routines.o  mps_grand.o
	$(COMPILER) $(FLAGS) eigensolver.o eigensolver_routines.o open_dmrg.o bose_hubbard.o mps_grand.o $(LIB) -o mps_grand.x

mps_grand.o: mps_grand.f90
	$(COMPILER)  -c mps_grand.f90

eigensolver.o: eigensolver.f
	$(COMPILER)  -c eigensolver.f

eigensolver_routines.o: eigensolver_routines.f90
	$(COMPILER)  -c eigensolver_routines.f90

open_dmrg.o: open_dmrg.f90
	$(COMPILER)  -c open_dmrg.f90

bose_hubbard.o: bose_hubbard.f90
	$(COMPILER)  -c bose_hubbard.f90

clean:
	rm -rf *.o *.mod *.dat mps_grand.x *~
