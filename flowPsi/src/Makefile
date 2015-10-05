include ../flowpsi.conf
include ../revision.conf
include $(LOCI_BASE)/Loci.conf

#COPTLESS = -O0 -g
#COPT=$(COPTLESS)


INCLUDES = -I../include


TARGET = flowpsi

OBJS =  \
	main.o gridTools.o gridSetup.o gridReader.o gridDynamic.o \
	gridAxisymmetry.o gridCoriolis.o gridComponents.o\
	bcInputs.o bcImpermeable.o bcInflow.o bcOutflow.o bcTools.o\
	bcIsentropic.o bcFixedMass.o bcViscousWall.o \
	solverTimestep.o solverNewton.o solverEoS.o solverIC.o \
	solverLinearSystem.o solverResidual.o solverRestart.o\
	solverTurbulence.o solverScalarTransport.o solverPeriodic.o \
	solverSetup.o solverGravity.o solverAveraging.o solverInterpolate.o \
	IOplot.o IOintegrated.o IOgrid.o IOprobe.o\
	linearSolverSGS.o linearSolverFSGS.o linearSolverLSGS.o \
	fluxInviscid.o fluxLowDiss.o fluxViscous.o



LOCAL_LIBS = $(LIBXML2_LIBS) -L$(LOCI_BASE)/lib 

JUNK = *~  core debug .schedule .schedule-*

default:
	$(MAKE) $(TARGET)

all:
	$(MAKE) $(TARGET)


$(TARGET): $(OBJS) 
	$(LD) -o  $(TARGET)  $(RPATH)$(FLOWPSI_RPATH) $(OBJS) $(LOCAL_LIBS) $(LIBS) $(LDFLAGS)


FRC : 

LOCI_FILES = $(wildcard *.loci)
LOCI_LPP_FILES = $(LOCI_FILES:.loci=.cc)

clean:
	rm -fr $(OBJS) $(TARGET) $(JUNK)

distclean: 
	rm -f $(DEPEND_FILES)
	rm -fr $(OBJS) $(TARGET) $(JUNK) $(LOCI_LPP_FILES)


DEPEND_FILES=$(subst .o,.d,$(OBJS))


#include automatically generated dependencies
-include $(DEPEND_FILES)
