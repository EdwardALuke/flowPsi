FLOWPSI_BASE?=..
include $(FLOWPSI_BASE)/flowpsi.conf


include $(LOCI_BASE)/Loci.conf

INCLUDES = -I$(FLOWPSI_BASE)/include 


# Put objects in the module here
OBJS = util.o par_util.o particle.o lagrangian_boundary.o lagrangian_statoutput.o particle_model.o particle_base.o particle_init.o integrate_particle.o load_particle_file.o flowPsi_interface.o

MODULE_NAME = Particle

LOCAL_LIBS = 

JUNK = *~  core ti_files ii_files rii_files

LIB_OBJS=$(OBJS:.o=_lo.o)

$(MODULE_NAME)_m.so: $(LIB_OBJS)
	$(SHARED_LD) $(SHARED_LD_FLAGS) $(MODULE_NAME)_m.so $(LIB_FLAGS) $(LIB_OBJS)

FRC : 

all: $(MODULE_NAME)_m.so 

clean:
	rm -fr $(OBJS) $(LIB_OBJS) $(MODULE_NAME)_m.so $(JUNK)

LOCI_FILES = $(wildcard *.loci)
LOCI_LPP_FILES = $(LOCI_FILES:.loci=.cc)

distclean: clean
	rm -fr $(DEPEND_FILES)
	rm -fr $(OBJS) $(LIB_OBJS) $(TARGET) $(JUNK) $(LOCI_LPP_FILES) $(MODULE_NAME)_m.so 

DEPEND_FILES=$(subst .o,.d,$(OBJS))


#include automatically generated dependencies
-include $(DEPEND_FILES)
