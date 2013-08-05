# Intel compilers on Unix
#CC = icc
#CFLAGS = -gcc-version=400 -O2 -mp
#LIBS = -lm -L. -lnum
#F77 = ifort
#FFLAGS = -gcc-version=400 -O2 -mp -nofor_main


CFLAGS = -ggdb
CC = gcc
LIBS = -lgsl -lgslcblas -lm -L. -lnum

#F77 = gfortran
#FFLAGS = -ggdb

all:	lib elekn jetkn flarekn pointkn

lib:	libn

libn:	nrutil.c save.c
	$(CC) $(CFLAGS) -c nrutil.c
	$(CC) $(CFLAGS) -c save.c
	rm -rf libnum.a
	ar rv libnum.a nrutil.o save.o
	ranlib libnum.a
	rm -rf nrutil.o save.o

elekn:	elekn.c ercph.c lib
	$(CC) $(CFLAGS) elekn.c model.c FS.c tridag.c ercph.c -o elekn $(LIBS)
	rm -rf elekn.o

pointkn:  pointkn.c ercph.c lib
	$(CC) $(CFLAGS) pointkn.c model.c ercph.c shell.c locate.c linear.c -o pointkn $(LIBS)
	rm -rf pointkn.o

jetkn:	jetkn.c ercph.c lib
	$(CC) $(CFLAGS) jetkn.c model.c ercph.c shell.c locate.c linear.c -o jetkn $(LIBS)
	rm -rf jetkn.o

flarekn:	flarekn.c lib
	$(CC) $(CFLAGS) flarekn.c model.c -o flarekn $(LIBS)

avekn:	avekn.c
	$(CC) $(CFLAGS) avekn.c shell.c locate.c linear.c -o avekn $(LIBS)
	rm -rf avekn.o

test:	test.c lib
	$(CC) $(CFLAGS) test.c -o test -lm $(LIBS)
	rm -rf test.o
