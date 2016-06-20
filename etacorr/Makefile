#

G77 = gfortran
CXX = g++
CC = gcc
LD = g++ 
LDFLAGS  = -L./ -L$(ROOTSYS)/lib

INCLUDES      = -I. -I$(shell root-config --incdir) 
ROOTCFLAGS   := $(shell root-config --cflags) 
ROOTLIBS     := $(shell root-config --libs)
ROOTGLIBS    := $(shell root-config --new --glibs)
ROOTGLIBS    += -lNew -lGui
CERNLIBS      = -L$(CERN)/pro/lib -lkernlib -lpacklib -lmathlib
LIBS          = $(ROOTLIBS) -lGui -lHtml -lGuiHtml $(SYSLIBS)
GLIBS         = $(ROOTGLIBS) $(SYSLIBS) 

CFLAGS = -g -O -I. -ITPC 
CXXFLAGS = -stdlib=libc++ -std=c++11 ${CFLAGS} $(INCLUDES) 

ifeq ($(MACHTYPE),x86_64)
	ARCHFLAG = -m64
else
	ARCHFLAG = -m32
endif

#DBG = -g -fno-inline -Wno-write-strings $(ARCHFLAG)
DBG = $(CXXFLAGS) -g -Wno-write-strings $(ARCHFLAG)

LDFLAGS = -Qunused-arguments -stdlib=libc++

all: etacorr

#----------------------------------------------------------------------------------
etacorr: etacorr.o
	${CXX} etacorr.o \
		$(LDFLAGS) $(ARCHFLAG) $(INCLUDES) $(LDLIBS) $(LIBS) -o etacorr

etacorr.o: etacorr.cxx
	${CXX} -c $(DBG) $< $(INCLUDES)
	
#fluct_main.o: fluct_class.h triggerdefs.h fluct_common.h
