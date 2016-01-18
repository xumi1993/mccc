CC = gcc
LIBS = -lm -lsac -lsacio -lfftw3
SACINC = -I/usr/local/sac/include
SACLIB = -L/usr/local/sac/lib 

all: mccc

mccc:  mccc.c sacio.c detrend.c
	${CC} -o $@ $^ ${SACINC} ${LIBS} ${SACLIB}

clean:
	rm mccc
