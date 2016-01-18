CC = gcc
LIBS = -lm -lsac -lsacio -L/usr/local/sac/lib -lfftw3
INCLUDE = -I/usr/local/sac/include

all: mccc

BIN = ./

mccc:  mccc.c sacio.c detrend.c
	${CC} -o $@ $^ ${INCLUDE} ${LIBS}

clean:
	rm *.o mccc
