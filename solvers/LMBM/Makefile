FF = gfortran -O2 -fPIC
CC = gcc

all: shared_lib

lmbmex.o: lmbmex.f
	$(FF) -c lmbmex.f

lmbm.o: lmbm.f
	$(FF) -c lmbm.f

lmsub.o: lmsub.f
	$(FF) -c lmsub.f

matcal.o: matcal.f
	$(FF) -c matcal.f

objfunc.o: objfunc.c
	$(CC) -c objfunc.c -o objfunc.o

shared_lib: lmbm.o lmbmex.o lmsub.o matcal.o objfunc.o
	gfortran -O3 -fPIC -shared objfunc.o lmbm.o lmbmex.o lmsub.o matcal.o -o liblmbm.so

clean:
	rm -f *.so *.o
