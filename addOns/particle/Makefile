# Makefile SETUP:
# Set FLOWPSI_BASE to the directory where flowPsi is installed
# Set MODULE_NAME to the name of your module
# Set FILES to list '.loci' files that will be compiled into your module, or
# use $(AUTOMATIC_FILES)

# What is the name of the module given to loadModule
MODULE_NAME = particleFPSI
# This will by default include all '.loci' files in the current directory
FILES = $(AUTOMATIC_FILES)

###########################################################################
# No changes required below this line
###########################################################################

AUTOMATIC_FILES = $(wildcard *.c *.cc *.cpp *.loci)
AUTOMATIC_OBJS = $(subst .c,.o,$(subst .cc,.o,$(subst .cpp,.o, $(subst .loci,.o,$(FILES)))))

OBJS := $(AUTOMATIC_OBJS)

include $(FLOWPSI_BASE)/flowpsi.conf
include $(LOCI_BASE)/Loci.conf

INCLUDES = -I$(FLOWPSI_BASE)/include -I./include

LOCAL_LIBS = 

MAKEFILE_DIRECTORY = $(shell pwd | sed -e s?/.*/??)

JUNK = *~  core 

LIB_OBJS=$(OBJS:.o=_lo.o)

all: $(MODULE_NAME)_m.so 

$(MODULE_NAME)_m.so: $(LIB_OBJS)
	$(SHARED_LD) $(SHARED_LD_FLAGS) $(MODULE_NAME)_m.so $(LIB_FLAGS) $(LIB_OBJS)

install_local: ../../lib/$(MODULE_NAME)_m.so

../../lib/$(MODULE_NAME)_m.so: $(MODULE_NAME)_m.so 
	cp $(MODULE_NAME)_m.so ../../lib/$(MODULE_NAME)_m.so

install: $(MODULE_NAME)_m.so 
	cp $(MODULE_NAME)_m.so $(INSTALL_PATH)/lib/$(MODULE_NAME)_m.so
ifneq ($(wildcard include/.),)
	@echo Installing Include Files
	mkdir -p $(INSTALL_PATH)/include/$(MAKEFILE_DIRECTORY)
	cp $(wildcard include/*.h include/*.lh include/*.hpp) $(INSTALL_PATH)/include/$(MAKEFILE_DIRECTORY)
endif


clean:
	rm -fr $(OBJS) $(LIB_OBJS) $(MODULE_NAME)_m.so $(JUNK)

LOCI_FILES = $(wildcard *.loci)
LOCI_LPP_FILES = $(LOCI_FILES:.loci=.$(LPP_I_SUFFIX))

distclean: clean
	rm -f $(DEPEND_FILES) $(LOCI_LPP_FILES)

DEPEND_FILES=$(subst .o,.d,$(OBJS))

#include automatically generated dependencies
ifeq ($(filter $(MAKECMDGOALS),clean distclean ),)
-include $(DEPEND_FILES)
endif
