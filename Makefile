include flowpsi.conf
include test.conf
include version.conf

FLOWPSI_BASE = $(shell pwd)

default: flowpsi tools turbulence guide addOns

all: default

INSTALL_PATH = $(FLOWPSI_INSTALL_PATH)/$(FLOWPSI_INSTALL_DIR)

install: FRC
	mkdir -p $(INSTALL_PATH)
	mkdir -p $(INSTALL_PATH)/lib
	mkdir -p $(INSTALL_PATH)/bin
	cp flowpsi.conf $(INSTALL_PATH)/flowpsi.conf
	@sed -e "s:^FLOWPSI_GIT_INFO.*:FLOWPSI_GIT_INFO = $(FLOWPSI_GIT_INFO):g" \
	    -e "s:^FLOWPSI_GIT_BRANCH.*:FLOWPSI_GIT_BRANCH = $(FLOWPSI_GIT_BRANCH):g" \
	                      revision.conf > \
	         $(INSTALL_PATH)/version.conf
	$(MAKE) -C src INSTALL_PATH="$(INSTALL_PATH)" FLOWPSI_BASE="$(FLOWPSI_BASE)" install
	$(MAKE) -C tools INSTALL_PATH="$(INSTALL_PATH)" FLOWPSI_BASE="$(FLOWPSI_BASE)" install
	$(MAKE) -C turbulence INSTALL_PATH="$(INSTALL_PATH)" FLOWPSI_BASE="$(FLOWPSI_BASE)" install
	$(MAKE) -C addOns INSTALL_PATH="$(INSTALL_PATH)" FLOWPSI_BASE="$(FLOWPSI_BASE)" install
	$(MAKE) -C guide INSTALL_PATH="$(INSTALL_PATH)" FLOWPSI_BASE="$(FLOWPSI_BASE)" install
	$(MAKE) -C examples INSTALL_PATH="$(INSTALL_PATH)" FLOWPSI_BASE="$(FLOWPSI_BASE)" install
	chmod -R a+rX $(INSTALL_PATH)


.PHONEY: FRC flowpsi tools test turbulence guide addOns install 

setup: FRC
	mkdir -p lib; true
	mkdir -p bin; true

flowpsi: setup
	$(MAKE) -C src FLOWPSI_BASE="$(FLOWPSI_BASE)" install_local

turbulence: setup
	$(MAKE) -C turbulence FLOWPSI_BASE="$(FLOWPSI_BASE)" install_local

addOns: setup
	$(MAKE) -C addOns FLOWPSI_BASE="$(FLOWPSI_BASE)" install_local

tools: setup
	$(MAKE) -C tools FLOWPSI_BASE="$(FLOWPSI_BASE)"  install_local

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
	$(MAKE) -C guide -k LOCI_BASE="$(LOCI_BASE)" FLOWPSI_BASE="$(FLOWPSI_BASE)" clean

distclean: FRC
	$(MAKE) -C quickTest -k LOCI_BASE="$(LOCI_BASE)" FLOWPSI_BASE="$(FLOWPSI_BASE)" clean
	$(MAKE) -C src -k LOCI_BASE="$(LOCI_BASE)" FLOWPSI_BASE="$(FLOWPSI_BASE)" distclean
	$(MAKE) -C tools -k LOCI_BASE="$(LOCI_BASE)" FLOWPSI_BASE="$(FLOWPSI_BASE)" distclean
	$(MAKE) -C turbulence -k  LOCI_BASE="$(LOCI_BASE)" FLOWPSI_BASE="$(FLOWPSI_BASE)" distclean
	$(MAKE) -C addOns -k  LOCI_BASE="$(LOCI_BASE)" FLOWPSI_BASE="$(FLOWPSI_BASE)" distclean
	$(MAKE) -C guide -k LOCI_BASE="$(LOCI_BASE)" FLOWPSI_BASE="$(FLOWPSI_BASE)" distclean
	rm -fr lib bin output debug *~ include/*~ version.conf flowpsi.conf


tarball:
	version_string=$$(sed -e "s:^FLOWPSI_GIT_INFO.*:FLOWPSI_GIT_INFO = $(FLOWPSI_GIT_INFO):g" \
	                      -e "s:^FLOWPSI_GIT_BRANCH.*:FLOWPSI_GIT_BRANCH = $(FLOWPSI_GIT_BRANCH):g" \
	                      revision.conf); \
	git archive --format=tgz --prefix=flowPsi-$(FLOWPSI_VERSION_NAME)/ \
	            --add-virtual-file=flowPsi-$(FLOWPSI_VERSION_NAME)/version.conf:"$$version_string" \
	            -o flowPsi-$(FLOWPSI_VERSION_NAME).tgz HEAD
