FLOWPSI_BASE ?= ..
include $(FLOWPSI_BASE)/flowpsi.conf
include $(FLOWPSI_BASE)/revision.conf
include $(LOCI_BASE)/Loci.conf

COPT = $(COPTLESS)

INCLUDES = -I$(FLOWPSI_BASE)/include 

OTHER_APPS = ndiff adpt

LOCAL_LIBS =

JUNK = *~ core ti_files ii_files rii_files 

all: $(OTHER_APPS)


FRC:


ndiff: ndiff.o
	$(LD) -o ndiff ndiff.o $(LOCAL_LIBS) $(LIBS) $(LDFLAGS)

adpt: adpt.o 
	$(LD) -o adpt adpt.o $(LOCAL_LIBS) $(LIBS) $(LDFLAGS)

clean:
	rm -fr $(OBJS) $(JUNK) make.depend make.depend.bak *.o $(OTHER_APPS)

distclean:
	rm -f $(DEPEND_FILES)
	rm -fr $(OBJS) $(TARGET) $(JUNK) make.depend make.depend.bak *.o $(OTHER_APPS)

DEPEND_FILES=$(subst .o,.d,$(OBJS)) $(subst .o,.d,$(EOSOBJS))


#include automatically generated dependencies
-include $(DEPEND_FILES)
