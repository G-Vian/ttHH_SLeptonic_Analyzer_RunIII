#----------------------------------------------------------------------------
# Description: Makefile to build analyzers with correctionlib
#----------------------------------------------------------------------------

ifndef ROOTSYS
$(error *** Please set up Root)
endif

name    := ttHHULanalyzer

# Sub-directories
srcdir	:= src
tmpdir	:= tmp
libdir	:= lib
incdir	:= include
libgsl  := /cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/gsl/2.2.1-omkpbe2/include/
liblhad := /cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/lhapdf/6.2.1-omkpbe3/include/
liblhadl := /cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/lhapdf/6.2.1-omkpbe3/lib/
libmem  := ./TTH/MEIntegratorStandalone/libs
libmemc := ./TTH/MEIntegratorStandalone/libs/proclib

$(shell mkdir -p tmp)
$(shell mkdir -p lib)

ifdef verbose
AT 	:=
else
AT	:= @
endif

#-----------------------------------------------------------------------
# sources and objects
#-----------------------------------------------------------------------
header  := $(incdir)/tnm.h
linkdef := $(incdir)/linkdef.h
cinthdr := $(srcdir)/dictionary.h
cintsrc := $(srcdir)/dictionary.cc

appsrcs	:= $(wildcard *.cc)
appobjects	:= $(addprefix $(tmpdir)/,$(appsrcs:.cc=.o))
applications := $(appsrcs:.cc=)

ccsrcs	:= $(filter-out $(cintsrc),$(wildcard $(srcdir)/*.cc))
sources	:= $(ccsrcs) $(cintsrc)
objects	:= $(subst $(srcdir)/,$(tmpdir)/,$(sources:.cc=.o))

#-----------------------------------------------------------------------
# Compilers
#-----------------------------------------------------------------------
CXX     := g++
LINK	:= g++
CINT	:= rootcint

#-----------------------------------------------------------------------
# Correctionlib flags
#-----------------------------------------------------------------------
CORRECTION_FLAGS := $(shell correction config --cflags --ldflags --rpath)

#-----------------------------------------------------------------------
# Include paths
#-----------------------------------------------------------------------
CPPFLAGS:= -I. -I$(incdir) -I$(libgsl) -I$(liblhad) -I$(srcdir) \
           $(shell root-config --cflags) $(cppflags) \
           $(CORRECTION_FLAGS)

#-----------------------------------------------------------------------
# Compiler flags
#-----------------------------------------------------------------------
CXXFLAGS:= -c -g -O2 -ansi -Wall -pipe -fPIC

#-----------------------------------------------------------------------
# Linker
#-----------------------------------------------------------------------
LD	:= $(LINK) -Wl,-rpath,$(ROOTSYS)/lib

OS	:= $(shell uname -s)
ifeq ($(OS),Darwin)
    LDSHARED	:= $(LD) -dynamiclib
    LDEXT       := .dylib
else
    LDSHARED	:= $(LD) -shared
    LDEXT       := .so
endif

LDFLAGS := -g

LIBS	:=  $(shell root-config --libs) \
           -L$(libdir) -L$(libmem) -L$(libmemc) -L$(liblhadl) \
           -lMinuit -lMathCore -lTMVA -lRooFit -lcuba -lopenloops -lLHAPDF \
           $(CORRECTION_FLAGS)

sharedlib := $(libdir)/libtnm$(LDEXT)

#-----------------------------------------------------------------------
# Rules
#-----------------------------------------------------------------------
all:	$(sharedlib) $(applications) 

bin:	$(applications)

lib:	$(sharedlib)

$(applications)	: %	: $(tmpdir)/%.o  $(sharedlib)
	@echo "---> Linking $@"
	$(AT)$(LD) $(LDFLAGS) $< $(LIBS) -ltnm -o $@

$(appobjects)	: $(tmpdir)/%.o	: %.cc
	@echo "---> Compiling application `basename $<`" 
	$(AT)$(CXX) $(CXXFLAGS) $(CPPFLAGS)  $< -o $@

$(sharedlib)	: $(objects)
	@echo "---> Linking `basename $@`"
	$(AT)$(LDSHARED) $(LDFLAGS) -fPIC $(objects) $(LIBS) -o $@

$(objects)	: $(tmpdir)/%.o	: $(srcdir)/%.cc
	@echo "---> Compiling `basename $<`" 
	$(AT)$(CXX) $(CXXFLAGS) $(CPPFLAGS)  $< -o $@

$(cintsrc)  : $(header) $(linkdef)
	@echo "---> Generating dictionary `basename $@`"
	$(AT)$(CINT) -f $@ -c -I. -Iinclude -I$(ROOTSYS)/include $+
	$(AT)mv $(srcdir)/*.pcm $(libdir)

clean:
	rm -rf $(tmpdir)/* $(libdir)/* $(srcdir)/dictionary* $(applications)
