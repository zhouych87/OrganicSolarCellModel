#LIB=-L/opt/intel/cmkl/9.0/lib/em64t -lmkl_lapack -lmkl_em64t -lguide -pthread
LIB=-L/opt/intel/cmkl/10.0.5.025/lib/em64t -lmkl_intel_lp64 -lmkl_core -lmkl_sequential
SOURCE=transferintegral.o

#main: oc-transferintegral bin clean

oct:	
	ifort  -FR -lowercase -O3 -c transferintegral.f90
	ifort -o calv $(SOURCE)  $(LIB) 

transferintegral.o : transferintegral.f90
	ifort  -FR -lowercase -O3 -c transferintegral.f90

clean:
	rm transferintegral
	rm *.o

bin:
	cp calv ~/bin
