all: MMStest gridMotion onlineStats actuatorDisk

MMStest: FRC
	$(MAKE) -C MMStest FLOWPSI_BASE="$(FLOWPSI_BASE)" all     

gridMotion: FRC
	$(MAKE) -C gridMotion FLOWPSI_BASE="$(FLOWPSI_BASE)" all     

onlineStats: FRC
	$(MAKE) -C onlineStats FLOWPSI_BASE="$(FLOWPSI_BASE)" all     

onlineStats: FRC
	$(MAKE) -C actuatorDisk FLOWPSI_BASE="$(FLOWPSI_BASE)" all     

clean: FRC
	$(MAKE) -C MMStest FLOWPSI_BASE="$(FLOWPSI_BASE)" clean     
	$(MAKE) -C gridMotion FLOWPSI_BASE="$(FLOWPSI_BASE)" clean     
	$(MAKE) -C onlineStats FLOWPSI_BASE="$(FLOWPSI_BASE)" clean     
	$(MAKE) -C actuatorDisk FLOWPSI_BASE="$(FLOWPSI_BASE)" clean     


distclean: FRC
	$(MAKE) -C MMStest FLOWPSI_BASE="$(FLOWPSI_BASE)" distclean
	$(MAKE) -C gridMotion FLOWPSI_BASE="$(FLOWPSI_BASE)" distclean
	$(MAKE) -C onlineStats FLOWPSI_BASE="$(FLOWPSI_BASE)" distclean
	$(MAKE) -C actuatorDisk FLOWPSI_BASE="$(FLOWPSI_BASE)" distclean

FRC: 

.PHONY: all MMStest gridMotion onlineStats actuatorDisk FRC
