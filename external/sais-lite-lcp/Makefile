# Makefile for suftest and test

# options
CC						= gcc
#CXX						= g++
#OUTPUT_OPTION	= -o $@
CFLAGS                          = -ffast-math -O9 -funroll-loops -DNDEBUG
#CFLAGS				= -O3 -fomit-frame-pointer -funroll-loops
#CXXFLAGS			= -O3 -fomit-frame-pointer
CPPFLAGS			= -Wall -DNDEBUG
#CPPFLAGS			= -Wall
LDFLAGS				= 
LDLIBS				= 
#TARGET_ARCH		=
INSTALL_DIR=../..


# targets
.PHONY: all
all: suftest libsaislcp

install: libsaislcp sais-lcp.h
	cp libsaislcp.a $(INSTALL_DIR)/lib
	cp sais-lcp.h $(INSTALL_DIR)/include 

libsaislcp: sais-lcp.o
	ar rcs libsaislcp.a sais-lcp.o

suftest: sais-lcp.o suftest.o
test:
	$(CC) -O -g -Wall test.c sais-lcp.c -o test
	./test
	$(RM) test test.exe

distclean: clean
clean:
	$(RM) suftest suftest.exe test test.exe sais.o suftest.o

# dependencies
sais.o suftest.o: sais-lcp.h Makefile

