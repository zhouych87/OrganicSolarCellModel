# ***********************************************************************
f90comp=gfortran -g -Wall -Wextra -Warray-temporaries -Wconversion -fimplicit-none -fbacktrace -ffree-line-length-0 -fcheck=all -ffpe-trap=invalid,zero,overflow,underflow -finit-real=nan 
f90comp=gfortran
# ***********************************************************************
split:  fileiov9.o morse3d.o readitp.o splitv9.o 
	$(f90comp) splitv9.o fileiov9.o morse3d.o readitp.o -o splitv9

fileiov9.o:  fileiov9.f90
	$(f90comp) -c fileiov9.f90
morse3d.o: morse3d.f90
	$(f90comp) -c morse3d.f90
readitp.o: readitp.f90
	$(f90comp) -c readitp.f90
splitv9.o:  splitv9.f90
	$(f90comp) -c splitv9.f90
# ***********************************************************************
clean:
	rm  *.o  *.exe lib.mod
# ***********************************************************************

