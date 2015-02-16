all: SA KOmega

SA: FRC
	cd SA; $(MAKE) FLOWPSI_BASE="$(FLOWPSI_BASE)" all; cd ..     
	cd SA; $(MAKE) FLOWPSI_BASE="$(FLOWPSI_BASE)" install; cd ..     

KOmega: FRC
	cd KOmega; $(MAKE) FLOWPSI_BASE="$(FLOWPSI_BASE)" all; cd ..     
	cd KOmega; $(MAKE) FLOWPSI_BASE="$(FLOWPSI_BASE)" install; cd ..     

clean: FRC
	cd SA; $(MAKE) FLOWPSI_BASE="$(FLOWPSI_BASE)" clean; cd ..
	cd KOmega; $(MAKE) FLOWPSI_BASE="$(FLOWPSI_BASE)" clean; cd ..

distclean: FRC
	cd SA; $(MAKE) FLOWPSI_BASE="$(FLOWPSI_BASE)" distclean; cd ..
	cd KOmega; $(MAKE) FLOWPSI_BASE="$(FLOWPSI_BASE)" distclean; cd ..

FRC: 
