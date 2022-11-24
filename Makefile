CC=gcc
LD=g++
CPP=g++

all:	tcount

tcount:	counting.o
	$(LD) -std=c++17 -O3 counting.o -o tcount -fopenmp

counting.o:	counting.c
	$(CC) -c counting.c -O3 -fopenmp

clean:
	rm -f *.o tcount *~
