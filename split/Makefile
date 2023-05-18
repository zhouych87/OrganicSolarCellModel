# ***********************************************************************
f90comp=gfortran -g -Wall -Wextra -Warray-temporaries -Wconversion -fimplicit-none -fbacktrace -ffree-line-length-0 -fcheck=all -ffpe-trap=invalid,zero,overflow,underflow -finit-real=nan 
f90comp=gfortran
# ***********************************************************************
split:  fileio.o morse3d.o readitp.o split.o 
	$(f90comp) split.o fileio.o morse3d.o readitp.o -o oscsplitv14

fileio.o:  fileio.f90
	$(f90comp) -c fileio.f90
morse3d.o: morse3d.f90
	$(f90comp) -c morse3d.f90
readitp.o: readitp.f90
	$(f90comp) -c readitp.f90
split.o:  split.f90
	$(f90comp) -c split.f90
# ***********************************************************************
clean:
	rm  *.o  *.exe lib.mod
# ***********************************************************************

