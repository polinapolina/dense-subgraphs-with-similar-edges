CFLAGS=-O4 -DBREAKPOINTS -c -Wall -fpic
BINDIR=bin
CC=gcc

pseudopar: ${BINDIR}/pseudopar
${BINDIR}/pseudopar:
	${CC} ${CFLAGS} src/1.0/pseudopar.c -o ${BINDIR}/pseudopar.o
	${CC} -shared -o ${BINDIR}/lib_pseudopar.so ${BINDIR}/pseudopar.o
clean:
	rm -f ${BINDIR}/*
