

CC=gcc
INSTALL=install -c -m 644
CFLAGS=-fPIC -ansi -pedantic -Wall -Wextra -I. -Dlint -Wno-unused-parameter -Wno-unused-variable


LIBS=-lm -lc

prefix=/usr/local
bindir=$(prefix)/bin
libdir=$(prefix)/lib
incdir=$(prefix)/include
srcdir=src
SOURCE := $(wildcard $(srcdir)/*.c)
INCS=$(srcdir)/*.h
DATA=binEphem.405 observatories.dat
OBJ=$(SOURCE:.c=.o)

.PHONY:	all

all: liborbfit
	@echo $(OBJ)

%.o: %.c $(INCS)
	$(CC) -c -o $@ $< $(CFLAGS) 

liborbfit: $(OBJ) 
	$(CC) -shared -o liborbfit.so $(OBJ)


install: all
	$(INSTALL) liborbfit.so $(libdir)
	$(INSTALL) $(INCS) $(incdir)
	$(INSTALL) data/observatories.dat $(libdir)
	uname -m | grep -q 64 && ${INSTALL} data/binEphem.405_64 ${libdir}/binEphem.405 || ${INSTALL} data/binEphem.405_64 ${libdir}/binEphem.405


clean:
	\rm $(OBJ)
	\rm liborbfit.so
