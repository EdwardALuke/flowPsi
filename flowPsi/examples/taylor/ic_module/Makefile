# Set FLOWPSI_BASE to point to the FlowPsi installation directory to compile
# this module
#
# To use a module, compile the module and then put the resulting shared object
# file either in the installed lib directory with chem or some other place that
# is contained in the LD_LIBRARY_PATH.  Once this is done, the module can be
# by placing 'loadModule: MODULE_NAME' into the top of the vars file
#
# Makefile SETUP:
# Set FLOWPSI_BASE to the directory where FLOWPSI is installed
# Set MODULE_NAME to the name of your module
# Set OBJS to list '.o' files that will be compiled into your module


FLOWPSI_BASE?=/home/lush/FlowPsiCode/flowPsi

MODULE_NAME=taylorGreenIC

# Put objects in the module here
OBJS =\
	taylorGreenIC.o

###########################################################################
# No changes required below this line
###########################################################################
include $(FLOWPSI_BASE)/flowpsi.conf
include $(LOCI_BASE)/Loci.conf

INCLUDES = -I$(FLOWPSI_BASE)/include 
#uncomment this for a debugging compile
#COPT=-O0 -g 


LOCAL_LIBS = 

JUNK = *~  core ti_files ii_files rii_files

LIB_OBJS=$(OBJS:.o=_lo.o)

all: $(MODULE_NAME)_m.so 

$(MODULE_NAME)_m.so: $(LIB_OBJS)
	$(SHARED_LD) $(SHARED_LD_FLAGS) $(MODULE_NAME)_m.so $(LIB_FLAGS) $(LIB_OBJS)
FRC : 

	
clean:
	rm -fr $(OBJS) $(LIB_OBJS) $(MODULE_NAME)_m.so $(JUNK)

LOCI_FILES = $(wildcard *.loci)
LOCI_LPP_FILES = $(LOCI_FILES:.loci=.cc)

distclean: 
	rm $(DEPEND_FILES)
	rm -fr $(OBJS) $(LIB_OBJS) $(MODULE_NAME)_m.so $(JUNK) $(LOCI_LPP_FILES)

DEPEND_FILES=$(subst .o,.d,$(OBJS))


#include automatically generated dependencies
-include $(DEPEND_FILES)
