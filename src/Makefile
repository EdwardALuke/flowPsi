# Makefile SETUP:
# Set FLOWPSI_BASE to the directory where flowPsi is installed
# Set MODULE_NAME to the name of your module
# Set FILES to list '.loci' files that will be compiled into your module, or
# use $(AUTOMATIC_FILES) to compile all '.loci' files in the present directory

# Name of application
TARGET = flowpsi
# This will by default inclue all ',loci' and other c++ and C files in the linking step
FILES = $(AUTOMATIC_FILES)
LOCAL_LIBS = 
LOCAL_INCLUDES = 


###########################################################################
# No changes required below this line
###########################################################################

AUTOMATIC_FILES = $(wildcard *.c *.cc *.cpp *.loci)
AUTOMATIC_OBJS = $(subst .c,.o,$(subst .cc,.o,$(subst .cpp,.o, $(subst .loci,.o,$(FILES)))))

OBJS := $(AUTOMATIC_OBJS)

include $(FLOWPSI_BASE)/flowpsi.conf
include $(FLOWPSI_BASE)/version.conf
include $(LOCI_BASE)/Loci.conf


INCLUDES = -I$(FLOWPSI_BASE)/include $(FLOWPSI_VERSION_FLAGS) $(LOCAL_INCLUDES)


JUNK = *~  core 

default:
	$(MAKE) $(TARGET)


install_local: $(TARGET)
	cp $(TARGET) $(FLOWPSI_BASE)/bin

install: $(TARGET)
	cp $(TARGET) $(INSTALL_PATH)/bin/$(TARGET)
	@echo Installing Include Files
	mkdir -p $(INSTALL_PATH)/include
	cp $(wildcard ../include/*.h ../include/*.lh) $(INSTALL_PATH)/include


$(TARGET): $(OBJS) 
	$(LD) -o  $(TARGET)  $(RPATH)$(FLOWPSI_RPATH) $(OBJS) -L$(LOCI_BASE)/lib -lfvmadaptfunc $(LOCAL_LIBS) $(LIBS) $(LDFLAGS)

FRC : 

LOCI_FILES = $(wildcard *.loci)
LOCI_LPP_FILES = $(LOCI_FILES:.loci=.$(LPP_I_SUFFIX))

clean:
	rm -fr $(OBJS) $(TARGET) $(JUNK)

distclean: 
	rm -f $(DEPEND_FILES)
	rm -fr $(OBJS) $(TARGET) $(JUNK) $(LOCI_LPP_FILES)


DEPEND_FILES=$(subst .o,.d,$(OBJS))


#include automatically generated dependencies
ifeq ($(filter $(MAKECMDGOALS),clean distclean ),)
-include $(DEPEND_FILES)
endif
