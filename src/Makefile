CC = gcc
CPP = g++
OPT= -O2
CFLAGS = -std=c99 -Wall -I/usr/include/mpich2
CPPFLAGS = -Wall
DEF= -DPERF_PROFILE -DPERF_ASSERT
AR = ar
ARTAG = rcv

libmpit	: profile.o timer.o mpiwrapper.o mpiwrapper_f.o estimator.o parser.o reporter.o
	$(AR) $(ARTAG) libmpit.a mpiwrapper.o mpiwrapper_f.o timer.o profile.o estimator.o parser.o reporter.o
	$(RM) -f *.o
reporter.o: reporter.cpp
	$(CPP) $(CPPFLAGS) $(DEF) -c reporter.cpp
parser.o: parser.c
	$(CC) $(CFLAGS) $(DEF) -c parser.c
estimator.o: estimator.c
	$(CC) $(CFLAGS) $(DEF) -c estimator.c
profile.o : profile.c	
	$(CC) $(CFLAGS) $(DEF) -c profile.c
mpiwrapper.o: mpiwrapper.c
	$(CC) $(CFLAGS) $(DEF) -c mpiwrapper.c
mpiwrapper_f.o: mpiwrapper_f.c
	$(CC) $(CFLAGS) $(DEF) -c mpiwrapper_f.c 
timer.o : timer.c	
	$(CC) $(CFLAGS) -c timer.c
clean :
	rm -f *.o	rm -f libmpit.a
