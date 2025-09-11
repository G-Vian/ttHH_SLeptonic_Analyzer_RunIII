#----------------------------------------------------------------------------
# Description: Makefile to build analyzers
# Created:     Wed Jul 20 18:33:34 2022 by mkanalyzer.py v2.0.2 15-Apr-2019
# Author:      Shakespeare's ghost
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

# Set this equal to the @ symbol to suppress display of instructions
# while make executes
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

# Construct list of sources to be compiled into applications
appsrcs	:= $(wildcard *.cc)
appobjects	:= $(addprefix $(tmpdir)/,$(appsrcs:.cc=.o))

# Construct list of applications
applications := $(appsrcs:.cc=)

# Construct list of sources to be compiled into shared library
ccsrcs	:= $(filter-out $(cintsrc),$(wildcard $(srcdir)/*.cc))
sources	:= $(ccsrcs) $(cintsrc)
objects	:= $(subst $(srcdir)/,$(tmpdir)/,$(sources:.cc=.o))

# Display list of applications to be built
#say	:= $(shell echo "appsrcs:     $(appsrcs)" >& 2)
#say	:= $(shell echo "appobjects:  $(appobjects)" >& 2)
#say	:= $(shell echo "objects:     $(objects)" >& 2)
#$(error bye!) 

#-----------------------------------------------------------------------
# 	Define which compilers and linkers to use
#-----------------------------------------------------------------------
# If clang++ exists use it, otherwise use g++
####COMPILER:= $(shell which clang++)
####ifneq ($(COMPILER),)
####CXX     := clang++
####LINK	:= clang++
####else
CXX     := g++
LINK	:= g++
####endif 
CINT	:= rootcint

#-----------------------------------------------------------------------
# 	Define paths to be searched for C++ header files (#include ....)
#-----------------------------------------------------------------------
CPPFLAGS:= -I. -I$(incdir) -I$(libgsl) -I$(liblhad) -I$(srcdir) -I$(HOME)/.local/lib/python3.6/site-packages/correctionlib/include $(shell root-config --cflags) $(cppflags)
#^ included the correction lib location ^ 

# 	Define compiler flags to be used
#	-c		perform compilation step only 
#	-g		include debug information in the executable file
#	-O2		optimize
#	-ansi	require strict adherance to C++ standard
#	-Wall	warn if source uses any non-standard C++
#	-pipe	communicate via different stages of compilation
#			using pipes rather than temporary files

CXXFLAGS:= -c -g -O2 -ansi -Wall -pipe -fPIC

#	C++ Linker
#   set path to ROOT libraries (Mac OS workaround)
LD	:= $(LINK) -Wl,-rpath,$(ROOTSYS)/lib

OS	:= $(shell uname -s)
ifeq ($(OS),Darwin)
    LDSHARED	:= $(LD) -dynamiclib
    LDEXT       := .dylib
else
    LDSHARED	:= $(LD) -shared
    LDEXT       := .so
endif

#	Linker flags

LDFLAGS := -g

# 	Libraries

#LIBS	:=  $(shell root-config --libs) -L$(libdir) -L$(libmem) -L$(libmemc) -L$(liblhadl) -lMinuit -lMathCore -lTMVA -lRooFit -lLHAPDF
LIBS	:=  $(shell root-config --libs) -L$(libdir) -L$(libmem) -L$(libmemc) -L$(liblhadl) -lMinuit -lMathCore -lTMVA -lRooFit -lcuba -lopenloops -lLHAPDF
####LIBS	:=  $(shell root-config --libs) -L$(libdir) -L$(libmem) -L$(libmemc) -L$(liblhadl) -lMinuit -lMathCore -lTMVA -lRooFit -lLHAPDF
sharedlib := $(libdir)/libtnm$(LDEXT)

#-----------------------------------------------------------------------
#	Rules
#	The structure of a rule is
#	target : source
#		command
#	The command makes a target from the source. 
#	$@ refers to the target
#	$< refers to the source
#-----------------------------------------------------------------------
all:	$(sharedlib) $(applications) 

bin:	$(applications)

lib:	$(sharedlib)

# Syntax:
# list of targets : target pattern : source pattern

# Make applications depend on shared libraries to force the latter
# to be built first

$(applications)	: %	: $(tmpdir)/%.o  $(sharedlib)
	@echo "---> Linking $@"
	$(AT)$(LD) $(LDFLAGS) $< $(LIBS) -ltnm -o $@

$(appobjects)	: $(tmpdir)/%.o	: %.cc
	@echo "---> Compiling application `basename $<`" 
	$(AT)$(CXX) $(CXXFLAGS) $(CPPFLAGS)  $< -o $@ # >& $*.FAILED
	@rm -rf $*.FAILED

$(sharedlib)	: $(objects)
	@echo "---> Linking `basename $@`"
	$(AT)$(LDSHARED) $(LDFLAGS) -fPIC $(objects) $(LIBS) -o $@

$(objects)	: $(tmpdir)/%.o	: $(srcdir)/%.cc
	@echo "---> Compiling `basename $<`" 
	$(AT)$(CXX) $(CXXFLAGS) $(CPPFLAGS)  $< -o $@ # >& $*.FAILED
	$(AT)rm -rf $*.FAILED

$(cintsrc)  : $(header) $(linkdef)
	@echo "---> Generating dictionary `basename $@`"
	$(AT)$(CINT) -f $@ -c -I. -Iinclude -I$(ROOTSYS)/include $+
	$(AT)mv $(srcdir)/*.pcm $(libdir)

# 	Define clean up rules
clean   :
	rm -rf $(tmpdir)/* $(libdir)/* $(srcdir)/dictionary* $(applications)
