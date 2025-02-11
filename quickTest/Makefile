
default: inviscid viscous restart
	cat inviscid/TestResults viscous/TestResults restart/TestResults > TestResults

all: default

.phony: FRC

inviscid: FRC
	$(MAKE) -C inviscid LOCI_BASE="$(LOCI_BASE)" FLOWPSI_BASE="$(FLOWPSI_BASE)" 

viscous: FRC
	$(MAKE) -C viscous LOCI_BASE="$(LOCI_BASE)" FLOWPSI_BASE="$(FLOWPSI_BASE)" 

restart: FRC
	$(MAKE) -C restart  LOCI_BASE="$(LOCI_BASE)" FLOWPSI_BASE="$(FLOWPSI_BASE)" 

clean: FRC
	$(MAKE) -C inviscid LOCI_BASE="$(LOCI_BASE)" FLOWPSI_BASE="$(FLOWPSI_BASE)" clean 
	$(MAKE) -C viscous LOCI_BASE="$(LOCI_BASE)" FLOWPSI_BASE="$(FLOWPSI_BASE)" clean ; cd ..
	$(MAKE) -C restart LOCI_BASE="$(LOCI_BASE)" FLOWPSI_BASE="$(FLOWPSI_BASE)" clean ; cd ..
	rm -f TestResults *~

FRC :

