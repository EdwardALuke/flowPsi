all: MMStest gridMotion onlineStats actuatorDisk particle
install_local: MMStest gridMotion onlineStats actuatorDisk particle

MMStest: FRC
	$(MAKE) -C MMStest FLOWPSI_BASE="$(FLOWPSI_BASE)" install_local     

gridMotion: FRC
	$(MAKE) -C gridMotion FLOWPSI_BASE="$(FLOWPSI_BASE)" install_local     

onlineStats: FRC
	$(MAKE) -C onlineStats FLOWPSI_BASE="$(FLOWPSI_BASE)" install_local

actuatorDisk: FRC
	$(MAKE) -C actuatorDisk FLOWPSI_BASE="$(FLOWPSI_BASE)" install_local

particle: FRC
	$(MAKE) -C particle FLOWPSI_BASE="$(FLOWPSI_BASE)" install_local

install: FRC
	$(MAKE) -C MMStest FLOWPSI_BASE="$(FLOWPSI_BASE)" INSTALL_PATH="$(INSTALL_PATH) install
	$(MAKE) -C gridMotion FLOWPSI_BASE="$(FLOWPSI_BASE)" INSTALL_PATH="$(INSTALL_PATH) install
	$(MAKE) -C onlineStats FLOWPSI_BASE="$(FLOWPSI_BASE)" INSTALL_PATH="$(INSTALL_PATH) install
	$(MAKE) -C actuatorDisk FLOWPSI_BASE="$(FLOWPSI_BASE)" INSTALL_PATH="$(INSTALL_PATH) install
	$(MAKE) -C particle FLOWPSI_BASE="$(FLOWPSI_BASE)" INSTALL_PATH="$(INSTALL_PATH) install

clean: FRC
	$(MAKE) -C MMStest FLOWPSI_BASE="$(FLOWPSI_BASE)" clean
	$(MAKE) -C gridMotion FLOWPSI_BASE="$(FLOWPSI_BASE)" clean
	$(MAKE) -C onlineStats FLOWPSI_BASE="$(FLOWPSI_BASE)" clean
	$(MAKE) -C actuatorDisk FLOWPSI_BASE="$(FLOWPSI_BASE)" clean
	$(MAKE) -C particle FLOWPSI_BASE="$(FLOWPSI_BASE)" clean


distclean: FRC
	$(MAKE) -C MMStest FLOWPSI_BASE="$(FLOWPSI_BASE)" distclean
	$(MAKE) -C gridMotion FLOWPSI_BASE="$(FLOWPSI_BASE)" distclean
	$(MAKE) -C onlineStats FLOWPSI_BASE="$(FLOWPSI_BASE)" distclean
	$(MAKE) -C actuatorDisk FLOWPSI_BASE="$(FLOWPSI_BASE)" distclean
	$(MAKE) -C particle FLOWPSI_BASE="$(FLOWPSI_BASE)" distclean

FRC: 

.PHONY: all MMStest gridMotion onlineStats actuatorDisk particle FRC
