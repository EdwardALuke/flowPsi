# Set FLOWPSI_BASE to point to the FLOWPSI installation directory to compile
# this module
#
# To use a module, compile the module and then put the resulting shared object
# file either in the installed lib directory with flowpsi or some other place that
# is contained in the LD_LIBRARY_PATH.  Once this is done, the module can be 
# accessed by placing 'loadModule: MODULE_NAME' into the top of the vars file
#
# Makefile SETUP:
# Set FLOWPSI_BASE to the directory where FLOWPSI is installed
# Set MODULE_NAME to the name of your module
# Set FILES to list '.loci' files that will be compiled into your module, or
# use $(AUTOMATIC_FILES) to compile all '.loci' files in the present directory

#Points to FLOWPSI installation directory
FLOWPSI_BASE ?= ../../

# What is the name of the module given to loadModule
MODULE_NAME = MMStest
# This will by default include all '.loci' files in the current directory
FILES = $(AUTOMATIC_FILES)


###########################################################################
# No changes required below this line
###########################################################################

AUTOMATIC_FILES = $(wildcard *.loci)
AUTOMATIC_OBJS = $(subst .c,.o,$(subst .cc,.o,$(subst .loci,.o,$(FILES))))

OBJS := $(AUTOMATIC_OBJS)

include $(FLOWPSI_BASE)/flowpsi.conf
include $(LOCI_BASE)/Loci.conf

INCLUDES = -I$(FLOWPSI_BASE)/include -I./include
#uncomment this for a debugging compile
#COPT=-O0 -g 


LOCAL_LIBS = 

JUNK = *~  core ti_files ii_files rii_files

LIB_OBJS=$(OBJS:.o=_lo.o)

all: $(MODULE_NAME)FPSI_m.so 

$(MODULE_NAME)FPSI_m.so: $(LIB_OBJS)
	$(SHARED_LD) $(SHARED_LD_FLAGS) $(MODULE_NAME)FPSI_m.so $(LIB_FLAGS) $(LIB_OBJS)
	mv $(MODULE_NAME)FPSI_m.so ../../lib/$(MODULE_NAME)FPSI_m.so

FRC : 

install: all
	cp $(MODULE_NAME)FPSI_m.so $(FLOWPSI_BASE)/lib/$(MODULE_NAME)FPSI_m.so

clean:
	rm -fr $(OBJS) $(LIB_OBJS) $(MODULE_NAME)FPSI_m.so $(JUNK)

LOCI_FILES = $(wildcard *.loci)
LOCI_LPP_FILES = $(LOCI_FILES:.loci=.cc)

distclean: 
	rm -fr $(DEPEND_FILES)
	rm -fr $(OBJS) $(LIB_OBJS) $(MODULE_NAME)FPSI_m.so $(JUNK) $(LOCI_LPP_FILES)

DEPEND_FILES=$(subst .o,.d,$(OBJS))

#include automatically generated dependencies
-include $(DEPEND_FILES)
