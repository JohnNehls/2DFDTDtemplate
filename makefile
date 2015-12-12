CC=g++
# CFLAGS=-Wall -I. -fopenmp -O3 -pedantic-errors -Wno-write-strings
CFLAGS=-Wall -I. -fopenmp -O3  -Wno-write-strings -std=c++11
LIBS=-lm

# SRCS=gridtmz.c materialtmz.c pmltmz.c printData.c sources.c updatetmz.c inputtmz.c main.c setup.c
# OBJS= $(SRCS:.c=.o)
#ALLDEFS  -DPOINT_SOURCE -DOPENMP -DDRUDE-DPLANE_SOURCE -DCENTER_RICKER_SOURCE -DLEFT_RICKER_SOURCE
# STDDEFS= -DPOINT_SOURCE -DOPENMP
# DRU = -DDRUDE

make:
	mpic++ -o mpi.exe $(CFLAGS) mpiCode.cpp

debug:
	mpic++ -o mpi.exe -DDEBUG $(CFLAGS) mpiCode.cpp

clean:
	rm *.out
	rm *.exe

tidy:
	rm *~
	rm *#
sample:
	mpic++ -o sample.exe $(CFLAGS) sampleProgram.cpp

all:
	rm *.out
	rm *.exe
	mpic++ -o mpi.exe $(CFLAGS) mpiCode.cpp




