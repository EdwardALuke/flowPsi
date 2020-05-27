include ../flowpsi.conf
include ../revision.conf
include $(LOCI_BASE)/Loci.conf

FVMADAPT_LIB = -lfvmadaptfunc

#COPTLESS = -O0 -g
#COPT=$(COPTLESS)


INCLUDES = -I../include


TARGET = flowpsi

OBJS =  \
	main.o gridTools.o gridSetup.o gridReader.o gridDynamic.o \
	gridAxisymmetry.o gridCoriolis.o gridComponents.o gridMRF.o\
	bcInputs.o bcTurboInterface.o bcTools.o bcInterface.o \
	bcPrescribed.o bcImpermeable.o bcInflow.o bcOutflow.o \
	bcIsentropic.o bcFixedMass.o bcFixedMassOutflow.o bcViscousWall.o \
	solverTimestep.o solverRungeKutta.o solverNewton.o solverEoS.o solverIC.o \
	solverLinearSystem.o solverResidual.o solverRestart.o\
	solverTurbulence.o solverScalarTransport.o solverPeriodic.o \
	solverSetup.o solverGravity.o solverAveraging.o solverInterpolate.o \
	solverICInterpolate.o \
	IOplot.o IOintegrated.o IOgrid.o IOprobe.o IOClipSurfaces.o \
	linearSolverSGS.o linearSolverFSGS.o linearSolverLSGS.o \
	fluxInviscid.o fluxLowDiss.o fluxViscous.o



LOCAL_LIBS = $(LIBXML2_LIBS) -L$(LOCI_BASE)/lib $(FVMADAPT_LIB)

JUNK = *~  core debug .schedule .schedule-*

default:
	$(MAKE) $(TARGET)

all:
	$(MAKE) $(TARGET)


$(TARGET): $(OBJS) 
	$(LD) -o  $(TARGET)  $(RPATH)$(FLOWPSI_RPATH) $(OBJS) $(LOCAL_LIBS) $(LIBS) $(LDFLAGS)
	mv flowpsi ../bin/flowpsi

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
