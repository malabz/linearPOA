CC?=gcc
CPP = g++

STATIC_LINK = 1

ifeq ($(STATIC_LINK), 1)
	STATIC_ARG = -static
else
	STATIC_ARG =
endif
WFLAGS= -Wall -Wextra  -DNDEBUG -Wno-ignored-qualifiers
#WFLAGS=-Wall -Werror -Wno-unused-function
OPT_FLAGS= -ffast-math -funroll-loops -m64 -D_FILE_OFFSET_BITS=64 -I../include/boost/include -I../include/threadpool/include #-fomit-frame-pointer
CFLAGS = $(WFLAGS) $(OPT_FLAGS) $(STATIC_ARG)
CPPFLAGS = ${CFLAGS} -std=c++17 $(OPT_FLAGS)
ifeq ($(STATIC_LINK), 1)
	LDFLAGS = -Wl,--whole-archive -lpthread -Wl,--no-whole-archive
else
	LDFLAGS = -pthread
endif
SYSTEM?=$(shell uname -s)
LOCATION = $(shell pwd)
MAKE_PLACE = $(LOCATION)/$(firstword $(MAKEFILE_LIST))
DIR = $(shell dirname $(MAKE_PLACE))

PREFIX = /usr/local
BINDIR = $(PREFIX)/bin

THREADS = 1

LIBOBJ_MAIN = dp.o hirschberg.o mtset.o mtxutl.o seqio.o poa.o linearpoa.o

TARGETS = linearpoa liblinearpoa.a
##

TARGETS_ALL = $(TARGETS)

DEBUG = 0
DNA = 0
TERMINATOR = 1
DEFINES = -DDEBUG=$(DEBUG)
DEFINES += -DTERMINATOR=$(TERMINATOR)
CFLAGS += $(DEFINES)
ifeq ($(DEBUG), 1)
	CFLAGS += -g -O0 -pg
else
	CFLAGS += -O3
endif
##

INPUT = input.txt

##

all:${TARGETS_ALL}

linearpoa: main.cpp $(LIBOBJ_MAIN)
	$(CPP) -o linearpoa main.cpp $(LIBOBJ_MAIN) $(CPPFLAGS) $(LDFLAGS) $(OPT_FLAGS)

liblinearpoa.a: $(LIBOBJ_MAIN)
	$(AR) rcs $@ $(LIBOBJ_MAIN)

clean:
	${RM} *.o ${TARGETS} *.a

install: all
	mkdir -p $(BINDIR)
	cp $(TARGETS) $(BINDIR)
	@echo Install to $(BINDIR) done.

###

%.o: %.cpp
	$(CPP) $(CPPFLAGS) $(OPT_FLAGS) -c -o $@ $<

%.o:%.c
	 $(CC) $(CFLAGS) -c -o $@ $<

%.64.o:%.c
	 $(CC) $(CFLAGS) -c -o $@ $<

##
