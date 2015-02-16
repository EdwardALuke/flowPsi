
default: inviscid viscous
	cat inviscid/TestResults viscous/TestResults > TestResults

all: default

inviscid: FRC
	cd inviscid && $(MAKE) LOCI_BASE="$(LOCI_BASE)" FLOWPSI_BASE="$(FLOWPSI_BASE)" ; cd ..

viscous: FRC
	cd viscous && $(MAKE) LOCI_BASE="$(LOCI_BASE)" FLOWPSI_BASE="$(FLOWPSI_BASE)" ; cd ..; 

clean: FRC
	cd inviscid && $(MAKE) LOCI_BASE="$(LOCI_BASE)" FLOWPSI_BASE="$(FLOWPSI_BASE)" clean ; cd ..
	cd viscous && $(MAKE) LOCI_BASE="$(LOCI_BASE)" FLOWPSI_BASE="$(FLOWPSI_BASE)" clean ; cd ..
	rm -f TestResults *~

FRC :

