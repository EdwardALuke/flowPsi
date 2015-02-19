include flowpsi.conf
include test.conf
include revision.conf

FLOWPSI_BASE = $(shell pwd)

default: flowpsi tools turbulence particle

all: default

.PHONEY: FRC flowpsi tools test turbulence particle

setup: FRC
	mkdir -p lib
	mkdir -p bin
	rm -f bin/flowpsi
	ln -s $(FLOWPSI_BASE)/src/flowpsi bin/flowpsi
	rm -f bin/ndiff
	ln -s $(FLOWPSI_BASE)/tools/ndiff bin/ndiff

flowpsi: FRC
	cd src; $(MAKE) LOCI_BASE="$(LOCI_BASE)" all

turbulence: setup
	cd turbulence; $(MAKE) FLOWPSI_BASE="$(FLOWPSI_BASE)" all

particle: setup
	cd particle; $(MAKE) FLOWPSI_BASE="$(FLOWPSI_BASE)" install

tools: FRC
	cd tools; $(MAKE) FLOWPSI_BASE="$(FLOWPSI_BASE)"  all

test: flowpsi tools turbulence setup
	cd quickTest; $(MAKE) LOCI_BASE="$(LOCI_BASE)" FLOWPSI_BASE="$(FLOWPSI_BASE)" all; 
	grep -i pass quickTest/TestResults && ( grep -i fail quickTest/TestResults && exit 1 || echo Test Success )

FRC : 

testclean: FRC
	cd quickTest; $(MAKE) LOCI_BASE="$(LOCI_BASE)" FLOWPSI_BASE="$(FLOWPSI_BASE)" clean; cd ..

clean: FRC
	cd quickTest; $(MAKE) LOCI_BASE="$(LOCI_BASE)" FLOWPSI_BASE="$(FLOWPSI_BASE)" clean; cd ..
	cd src; $(MAKE) LOCI_BASE="$(LOCI_BASE)" FLOWPSI_BASE="$(FLOWPSI_BASE)" clean; cd ..	
	cd tools; $(MAKE) LOCI_BASE="$(LOCI_BASE)" FLOWPSI_BASE="$(FLOWPSI_BASE)" clean; cd ..
	cd turbulence; $(MAKE) LOCI_BASE="$(LOCI_BASE)" FLOWPSI_BASE="$(FLOWPSI_BASE)" clean; cd ..
	cd particle; $(MAKE) LOCI_BASE="$(LOCI_BASE)" FLOWPSI_BASE="$(FLOWPSI_BASE)" clean; cd ..

distclean: FRC
	cd quickTest; $(MAKE) LOCI_BASE="$(LOCI_BASE)" FLOWPSI_BASE="$(FLOWPSI_BASE)" clean; cd ..
	cd src; $(MAKE) LOCI_BASE="$(LOCI_BASE)" FLOWPSI_BASE="$(FLOWPSI_BASE)" distclean; cd ..
	cd tools; $(MAKE) LOCI_BASE="$(LOCI_BASE)" FLOWPSI_BASE="$(FLOWPSI_BASE)" distclean; cd ..
	cd turbulence; $(MAKE) LOCI_BASE="$(LOCI_BASE)" FLOWPSI_BASE="$(FLOWPSI_BASE)" distclean; cd ..
	cd particle; $(MAKE) LOCI_BASE="$(LOCI_BASE)" FLOWPSI_BASE="$(FLOWPSI_BASE)" distclean; cd ..
	rm -fr lib bin output debug *~ include/*~