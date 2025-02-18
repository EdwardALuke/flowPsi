# Set FLOWPSI_BASE to point to the FLOWPSI installation directory to compile
# this module

# To use a module, compile the module and then put the resulting shared object
# file either in the installed lib directory with flowpsi or some other place 
# that is contained in the LD_LIBRARY_PATH.  Once this is done, the module 
# can be loaded into flowpsi using the -load_module option.  For example
# chem -load_module extensions case_name
# where extensions above is the MODULE_NAME given below
# multiple modules can be loaded with multiple -load_module options
include $(FLOWPSI_BASE)/flowpsi.conf
include $(LOCI_BASE)/Loci.conf
include $(FLOWPSI_BASE)/version.conf

INCLUDES = -DUSE_LOCI_TYPES -I$(FLOWPSI_BASE)/include -I./include $(FLOWPSI_VERSION_FLAGS)

# Put objects in the module here
OBJS =	approx.o \
	constraints.o \
	idw.o \
	nodeData.o \
	treeApprox.o \
	rotor.o \
	geomRep.o \
	couplingLevel.o \
	constrainedSurfaces.o \
	flowpsi_interface.o 

# List your object files here
VOGDEFOBJS =  vogdeform.o \
	approx.o \
	constraints.o \
	idw.o \
	nodeData.o \
	treeApprox.o \
	rotor.o \
	geomRep.o \
	couplingLevel.o \
	constrainedSurfaces.o

MODULE_NAME = gridMotionFPSI

LOCAL_LIBS = 

MAKEFILE_DIRECTORY = $(shell pwd | sed -e s?/.*/??)

JUNK = *~  core 

LIB_OBJS=$(OBJS:.o=_lo.o)

$(MODULE_NAME)_m.so: $(LIB_OBJS)
	$(SHARED_LD) $(SHARED_LD_FLAGS) $(MODULE_NAME)_m.so $(LIB_FLAGS) $(LIB_OBJS) $(LIBS)

install_local: ../../lib/$(MODULE_NAME)_m.so ../../bin/vogdeform

../../lib/$(MODULE_NAME)_m.so: $(MODULE_NAME)_m.so 
	cp $(MODULE_NAME)_m.so ../../lib/$(MODULE_NAME)_m.so

../../bin/vogdeform: vogdeform
	cp vogdeform ../../bin/vogdeform

all: $(MODULE_NAME)_m.so vogdeform

install: $(MODULE_NAME)_m.so 
	cp $(MODULE_NAME)_m.so $(INSTALL_PATH)/lib/$(MODULE_NAME)_m.so
	cp vogdeform $(INSTALL_PATH)/bin/vogdeform
ifeq ($(shell test -d include),0)
	@echo Installing Include Files
	mkdir -p $(INSTALL_PATH)/include/$(MAKEFILE_DIRECTORY)
	cp $(wildcard include/*.h include/*.lh include/*.hpp) $(INSTALL_PATH)/include/$(MAKEFILE_DIRECTORY)
endif

vogdeform: $(VOGDEFOBJS)
	$(LD) -o vogdeform $(VOGDEFOBJS) $(LOCAL_LIBS) $(LIBS) $(LDFLAGS)

FRC : 

clean:
	rm -fr $(OBJS) $(LIB_OBJS) $(MODULE_NAME)_m.so $(VOGDEFOBJS) $(JUNK) vogdeform

LOCI_FILES = $(wildcard *.loci)
LOCI_LPP_FILES = $(LOCI_FILES:.loci=.$(LPP_I_SUFFIX))

distclean: 
	rm -fr $(DEPEND_FILES)
	rm -fr $(OBJS) $(LIB_OBJS) $(MODULE_NAME)_m.so $(VOGDEFOBJS) $(JUNK) $(LOCI_LPP_FILES)

DEPEND_FILES=$(subst .o,.d,$(OBJS)) 


#include automatically generated dependencies
ifeq ($(filter $(MAKECMDGOALS),clean distclean ),)
-include $(DEPEND_FILES)
endif
