include flowpsi.conf
include test.conf
include revision.conf

FLOWPSI_BASE = $(shell pwd)

default: flowpsi tools turbulence particle guide addOns

all: default

install: all
	FLOWPSI_INSTALL_DIR=$(FLOWPSI_INSTALL_DIR) FLOWPSI_INSTALL_PATH=$(FLOWPSI_INSTALL_PATH) bash Install.bash

.PHONEY: FRC flowpsi tools test turbulence particle guide addOns install

setup: FRC
	mkdir -p lib
	mkdir -p bin

flowpsi: FRC
	$(MAKE) -C src LOCI_BASE="$(LOCI_BASE)" all

turbulence: setup
	$(MAKE) -C turbulence FLOWPSI_BASE="$(FLOWPSI_BASE)" all

addOns: setup
	$(MAKE) -C addOns FLOWPSI_BASE="$(FLOWPSI_BASE)" all

particle: setup
	$(MAKE) -C particle FLOWPSI_BASE="$(FLOWPSI_BASE)" install

tools: FRC
	$(MAKE) -C tools FLOWPSI_BASE="$(FLOWPSI_BASE)"  all

docs: FRC
	$(MAKE) -C guide FLOWPSI_BASE="$(FLOWPSI_BASE)" guide.pdf

test: flowpsi tools turbulence setup
	$(MAKE) -C quickTest LOCI_BASE="$(LOCI_BASE)" FLOWPSI_BASE="$(FLOWPSI_BASE)" all; 
	grep -i pass quickTest/TestResults && ( grep -i fail quickTest/TestResults && exit 1 || echo Test Success )

FRC : 

testclean: FRC
	$(MAKE) -C quickTest LOCI_BASE="$(LOCI_BASE)" FLOWPSI_BASE="$(FLOWPSI_BASE)" clean; cd ..

clean: FRC
	$(MAKE) -C quickTest -k LOCI_BASE="$(LOCI_BASE)" FLOWPSI_BASE="$(FLOWPSI_BASE)" clean
	$(MAKE) -C src -k LOCI_BASE="$(LOCI_BASE)" FLOWPSI_BASE="$(FLOWPSI_BASE)" clean
	$(MAKE) -C tools -k LOCI_BASE="$(LOCI_BASE)" FLOWPSI_BASE="$(FLOWPSI_BASE)" clean
	$(MAKE) -C turbulence -k  LOCI_BASE="$(LOCI_BASE)" FLOWPSI_BASE="$(FLOWPSI_BASE)" clean
	$(MAKE) -C addOns -k LOCI_BASE="$(LOCI_BASE)" FLOWPSI_BASE="$(FLOWPSI_BASE)" clean
	$(MAKE) -C particle -k LOCI_BASE="$(LOCI_BASE)" FLOWPSI_BASE="$(FLOWPSI_BASE)" clean
	$(MAKE) -C guide -k LOCI_BASE="$(LOCI_BASE)" FLOWPSI_BASE="$(FLOWPSI_BASE)" clean

distclean: FRC
	$(MAKE) -C quickTest -k LOCI_BASE="$(LOCI_BASE)" FLOWPSI_BASE="$(FLOWPSI_BASE)" clean
	$(MAKE) -C src -k LOCI_BASE="$(LOCI_BASE)" FLOWPSI_BASE="$(FLOWPSI_BASE)" distclean
	$(MAKE) -C tools -k LOCI_BASE="$(LOCI_BASE)" FLOWPSI_BASE="$(FLOWPSI_BASE)" distclean
	$(MAKE) -C turbulence -k  LOCI_BASE="$(LOCI_BASE)" FLOWPSI_BASE="$(FLOWPSI_BASE)" distclean
	$(MAKE) -C addOns -k  LOCI_BASE="$(LOCI_BASE)" FLOWPSI_BASE="$(FLOWPSI_BASE)" distclean
	$(MAKE) -C particle -k LOCI_BASE="$(LOCI_BASE)" FLOWPSI_BASE="$(FLOWPSI_BASE)" distclean
	$(MAKE) -C guide -k LOCI_BASE="$(LOCI_BASE)" FLOWPSI_BASE="$(FLOWPSI_BASE)" distclean
	rm -fr lib bin output debug *~ include/*~
