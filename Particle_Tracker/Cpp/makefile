#!/bin/tcsh
#
# Makefile for PARTRAC
#
CC = g++
CFLAGS = -Wall -fbounds-checking
FFLAGS =
OBJECTS = main.o motion.o init.o weighting.o
EXEC = partrac
OPTFLAGS = -O3

all : ${EXEC}

.PHONY : clean

clean::
	rm -f ${OBJECTS}
	rm -f ${EXEC}

${EXEC} : ${OBJECTS}
	${CC} ${CFLAGS} ${OPTFLAGS} -o ${EXEC} ${OBJECTS}

%.o : %.cpp
	$(CC) $(CF) $(INCPATH) -c $<
