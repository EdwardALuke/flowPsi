include flowpsi.conf
include test.conf
include version.conf

FLOWPSI_BASE = $(shell pwd)

default: flowpsi tools turbulence guide addOns particle 

all: default

install: all
	FLOWPSI_INSTALL_DIR=$(FLOWPSI_INSTALL_DIR) FLOWPSI_INSTALL_PATH=$(FLOWPSI_INSTALL_PATH) bash Install.bash

.PHONEY: FRC flowpsi tools test turbulence guide addOns install particle 

setup: FRC
	mkdir -p lib; true
	mkdir -p bin; true

flowpsi: setup
	$(MAKE) -C src LOCI_BASE="$(LOCI_BASE)" all

turbulence: setup
	$(MAKE) -C turbulence FLOWPSI_BASE="$(FLOWPSI_BASE)" all

addOns: setup
	$(MAKE) -C addOns FLOWPSI_BASE="$(FLOWPSI_BASE)" all

particle: setup
	$(MAKE) -C particle FLOWPSI_BASE="$(FLOWPSI_BASE)" install

tools: setup
	$(MAKE) -C tools FLOWPSI_BASE="$(FLOWPSI_BASE)"  all

docs: FRC
	$(MAKE) -C guide FLOWPSI_BASE="$(FLOWPSI_BASE)" flowPsiGuide.pdf

test: flowpsi tools turbulence
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
	rm -fr lib bin output debug *~ include/*~ version.conf flowpsi.conf


tarball:
	version_string=$$(sed -e "s:^FLOWPSI_GIT_INFO.*:FLOWPSI_GIT_INFO = $(FLOWPSI_GIT_INFO):g" \
	                      -e "s:^FLOWPSI_GIT_BRANCH.*:FLOWPSI_GIT_BRANCH = $(FLOWPSI_GIT_BRANCH):g" \
	                      revision.conf); \
	git archive --format=tgz --prefix=flowPsi-$(FLOWPSI_VERSION_NAME)/ \
	            --add-virtual-file=flowPsi-$(FLOWPSI_VERSION_NAME)/version.conf:"$$version_string" \
	            -o flowPsi-$(FLOWPSI_VERSION_NAME).tgz HEAD
