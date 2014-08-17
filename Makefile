# ============================================================================
# Name        : Makefile
# Author      : amir
# Version     :
# Copyright   : Your copyright notice
# Description : Makefile for Hello World in Fortran
# ============================================================================

.PHONY: all clean run

# Change this line if you are using a different Fortran compiler
FORTRAN_COMPILER = gfortran
exec_file=    bin\run.exe
module_files= driver/linearsolvers.f90 driver/femgeometry.f90 driver/materialbehavior.f90 driver/femlibs.f90 
driver_file = driver/driver_fem_1d_fortran_eclipse.f90
compiled_module_files=fem_geometry.mod fem_libs.mod linearsolvers.mod material_behavior.mod
output_files=*.mod *.exe *.csv *.out *.txt fort.*
library_files=c:\lapack\liblapack.a c:\lapack\librefblas.a c:\lapack\libtmglib.a 


all: 
	del $(output_files) $(exec_file)
	$(FORTRAN_COMPILER) -O2 -g -o \
	$(exec_file) \
	$(module_files) $(driver_file)	$(library_files)
	

clean:
	del $(output_files) $(exec_file) ;
	
	
run: 
	del $(output_files) $(exec_file)
	$(FORTRAN_COMPILER) -O2 -g -o \
	$(exec_file) \
	$(module_files) $(driver_file) $(library_files)
	$(exec_file)
