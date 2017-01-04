all: MMStest gridMotion onlineStats 

MMStest: FRC
	$(MAKE) -C MMStest FLOWPSI_BASE="$(FLOWPSI_BASE)" all     

gridMotion: FRC
	$(MAKE) -C gridMotion FLOWPSI_BASE="$(FLOWPSI_BASE)" all     

onlineStats: FRC
	$(MAKE) -C onlineStats FLOWPSI_BASE="$(FLOWPSI_BASE)" all     

clean: FRC
	$(MAKE) -C MMStest FLOWPSI_BASE="$(FLOWPSI_BASE)" clean     
	$(MAKE) -C gridMotion FLOWPSI_BASE="$(FLOWPSI_BASE)" clean     
	$(MAKE) -C onlineStats FLOWPSI_BASE="$(FLOWPSI_BASE)" clean     


distclean: FRC
	$(MAKE) -C MMStest FLOWPSI_BASE="$(FLOWPSI_BASE)" distclean
	$(MAKE) -C gridMotion FLOWPSI_BASE="$(FLOWPSI_BASE)" distclean
	$(MAKE) -C onlineStats FLOWPSI_BASE="$(FLOWPSI_BASE)" distclean

FRC: 

.PHONY: all MMStest gridMotion onlineStats FRC
