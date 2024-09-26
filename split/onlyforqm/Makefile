# ***********************************************************************
f90comp=gfortran -g -Wall -Wextra -Warray-temporaries -Wconversion -fimplicit-none -fbacktrace -ffree-line-length-0 -fcheck=all -ffpe-trap=invalid,zero,overflow,underflow -finit-real=nan 
f90comp=gfortran
# ***********************************************************************
split:  fileio.o  split.o 
	$(f90comp) split.o fileio.o -o oscsplit

fileio.o:  fileio.f90
	$(f90comp) -c fileio.f90
split.o:  split.f90
	$(f90comp) -c split.f90
# ***********************************************************************
clean:
	rm  *.o  *.exe lib.mod
# ***********************************************************************

