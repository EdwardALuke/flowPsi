all: SA KOmega

SA: FRC
	$(MAKE) -C SA FLOWPSI_BASE="$(FLOWPSI_BASE)" all
	$(MAKE) -C FLOWPSI_BASE="$(FLOWPSI_BASE)" install

KOmega: FRC
	$(MAKE) -C KOmega FLOWPSI_BASE="$(FLOWPSI_BASE)" all
	$(MAKE) -C KOmega FLOWPSI_BASE="$(FLOWPSI_BASE)" install

clean: FRC
	$(MAKE) -C SA -k FLOWPSI_BASE="$(FLOWPSI_BASE)" clean
	$(MAKE) -C KOmega -k FLOWPSI_BASE="$(FLOWPSI_BASE)" clean

distclean: FRC
	$(MAKE) -C SA -k FLOWPSI_BASE="$(FLOWPSI_BASE)" distclean
	$(MAKE) -C KOmega -k FLOWPSI_BASE="$(FLOWPSI_BASE)" distclean

FRC: 

.PHONY: all SA KOmega clean distclean FRC
