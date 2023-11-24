# ***********************************************************************
#f90comp=gfortran -g -Wall -Wextra -Warray-temporaries -Wconversion -fimplicit-none -fbacktrace -ffree-line-length-0 -fcheck=all -ffpe-trap=invalid,zero,overflow,underflow -finit-real=nan 
#f90comp=gfortran
f90comp=mpiifort
#f90comp=gfortran -g -fcheck=all -Wall -fbacktrace
# ***********************************************************************
desgen:  fileio.o morse3d.o readitp.o desgen.o 
	$(f90comp) desgen.o fileio.o morse3d.o readitp.o -o desfftmpi

fileio.o:  fileio.f90
	$(f90comp) -c fileio.f90
morse3d.o: morse3d.f90
	$(f90comp) -c morse3d.f90
readitp.o: readitp.f90
	$(f90comp) -c readitp.f90
desgen.o:  desgen.f90
	$(f90comp) -c desgen.f90
# ***********************************************************************
clean:
	rm  *.o  *.exe lib.mod
# ***********************************************************************

