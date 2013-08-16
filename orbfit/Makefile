

CC=gcc

CFLAGS=-fPIC -ansi -pedantic -Wall -Wextra -I. -Dlint -Wno-unused-parameter -Wno-unused-variable

LIBS=-lm -lc

OBJ=covsrt.o \
	dms.o	\
	ephem_earth.o \
	gasdev.o \
	gaussj.o \
	lubksb.o \
	ludcmp.o \
	mrqcof_orbit.o \
	mrqmin_orbit.o \
	nrutil.o \
	orbfit1.o \
	orbfit2.o \
	ran1.o \
	transforms.o \
	orbfitmodule.o


DEPS=ephem_read.h ephem_types.h nrutil.h orbfit.h 

%o: %c $(DEPS) 
	$(CC) -c -o $@ $< $(CFLAGS) 

liborbfit: $(OBJ)
	$(CC) -shared -o liborbfit.so $(OBJ)

