CC = gcc
CPP = g++
OPT= -O2
CFLAGS = -std=c99 $(CPPFLAGS)
CPPFLAGS = $(OPT) -g -Wall -I/usr/include 
DEF= -DPERF_ASSERT -DALLTRACE -D_XOPEN_SOURCE=500
AR = ar
ARTAG = rcv

OBJS= timer.o mpiwrapper.o mpiwrapper_f.o estimator.o parser.o mpiname.o reporter.o

ALL: main.out packer.out analist.out
analist.out: analist.cpp
	$(CPP) $(CPPFLAGS) analist.cpp -o analist.out
packer.out: packer.cpp
	$(CPP) $(CPPFLAGS) -fopenmp $(DEF) packer.cpp -o packer.out
main.out: main.o libmpit
	g++ main.o -L. -L/usr/bin -lmpit -lmpich -o main.out
libmpit	: $(OBJS)
	$(AR) $(ARTAG) libmpit.a $(OBJS)
reporter.o: reporter.cpp
	$(CPP) $(CPPFLAGS) $(DEF) -c reporter.cpp
parser.o: parser.c
	$(CC) $(CFLAGS) $(DEF) -c parser.c
estimator.o: estimator.c
	$(CC) $(CFLAGS) $(DEF) -c estimator.c
mpiwrapper.o: mpiwrapper.c
	$(CC) $(CFLAGS) $(DEF) -c mpiwrapper.c
mpiwrapper_f.o: mpiwrapper_f.c
	$(CC) $(CFLAGS) $(DEF) -c mpiwrapper_f.c 
timer.o : timer.c	
	$(CC) $(CFLAGS) -c timer.c
mpiname.o: mpiname.c mpiname.h
	$(CC) $(CFLAGS) -c mpiname.c
main.o: main.c
	$(CC) $(CFLAGS) -c main.c
clean :
	rm -f *.o	
	rm -f libmpit.a
	rm -f *.out
	rm -f *.log
	rm -f mpi_record.*
	rm -f mpipa.report.*
	rm -f mpipa.log.*
