# Makefile SETUP:
# Set FLOWPSI_BASE to the directory where flowPsi is installed
# Set APPS to list of applciations

APPS = ndiff adpt


include $(FLOWPSI_BASE)/flowpsi.conf
include $(FLOWPSI_BASE)/version.conf
include $(LOCI_BASE)/Loci.conf


INCLUDES = -I$(FLOWPSI_BASE)/include $(FLOWPSI_VERSION_FLAGS)

LOCAL_LIBS =

JUNK = *~ core 

LOCAL_APPS = $(addprefix $(FLOWPSI_BASE)/bin/,$(APPS))

install_local: $(LOCAL_APPS) 

INSTALL_APPS = $(addprefix $(INSTALL_PATH)/bin/,$(APPS))

install: $(INSTALL_APPS)

$(FLOWPSI_BASE)/bin/ndiff: ndiff
	cp ndiff $(FLOWPSI_BASE)/bin/ndiff

$(FLOWPSI_BASE)/bin/adpt: adpt
	cp adpt $(FLOWPSI_BASE)/bin

$(INSTALL_PATH)/bin/ndiff: ndiff
	cp ndiff $(INSTALL_PATH)/bin/ndiff

$(INSTALL_PATH)/bin/adpt: adpt
	cp adpt $(INSTALL_PATH)/bin/adpt

FRC:


ndiff: ndiff.o
	$(LD) -o ndiff ndiff.o $(LOCAL_LIBS) $(LIBS) $(LDFLAGS)

adpt: adpt.o 
	$(LD) -o adpt adpt.o $(LOCAL_LIBS) $(LIBS) $(LDFLAGS)


clean:
	rm -fr $(OBJS) $(JUNK) $(APPS)

distclean: clean
	rm -f $(DEPEND_FILES)

DEPEND_FILES=$(subst .o,.d,$(OBJS)) $(subst .o,.d,$(EOSOBJS))


#include automatically generated dependencies
ifeq ($(filter $(MAKECMDGOALS),clean distclean ),)
-include $(DEPEND_FILES)
endif
