all: SA KOmega DHRL

SA: FRC
	$(MAKE) -C SA FLOWPSI_BASE="$(FLOWPSI_BASE)" all
	$(MAKE) -C SA FLOWPSI_BASE="$(FLOWPSI_BASE)" install

KOmega: FRC
	$(MAKE) -C KOmega FLOWPSI_BASE="$(FLOWPSI_BASE)" all
	$(MAKE) -C KOmega FLOWPSI_BASE="$(FLOWPSI_BASE)" install

DHRL: FRC
	$(MAKE) -C DHRL FLOWPSI_BASE="$(FLOWPSI_BASE)" all
	$(MAKE) -C DHRL FLOWPSI_BASE="$(FLOWPSI_BASE)" install

clean: FRC
	$(MAKE) -C SA -k FLOWPSI_BASE="$(FLOWPSI_BASE)" clean
	$(MAKE) -C KOmega -k FLOWPSI_BASE="$(FLOWPSI_BASE)" clean
	$(MAKE) -C DHRL -k FLOWPSI_BASE="$(FLOWPSI_BASE)" clean

distclean: FRC
	$(MAKE) -C SA -k FLOWPSI_BASE="$(FLOWPSI_BASE)" distclean
	$(MAKE) -C KOmega -k FLOWPSI_BASE="$(FLOWPSI_BASE)" distclean
	$(MAKE) -C DHRL -k FLOWPSI_BASE="$(FLOWPSI_BASE)" distclean

FRC: 

.PHONY: all SA KOmega DHRL clean distclean FRC
