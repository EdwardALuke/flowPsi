
MODULE_CASES=$(filter-out moduleGuide/. scalarTransport/., $(wildcard */.))

.phony: FRC setup copy

setup: FRC
	mkdir -p $(INSTALL_PATH)/examples

copy :$(MODULE_CASES) scalarTransport moduleGuide setup

install: copy 


scalarTransport: FRC
	mkdir -p $(INSTALL_PATH)/examples/moduleGuide/scalarTransport
	mkdir -p $(INSTALL_PATH)/examples/moduleGuide/scalarTransport/MTR
	mkdir -p $(INSTALL_PATH)/examples/moduleGuide/scalarTransport/tracers
	cp $(wildcard scalarTransport/MTR/*.vars scalarTransport/MTR/*.vog \
	              scalarTransport/MTR/*.loci scalarTransport/MTR/*.lh \
	              scalarTransport/MTR/*.h) Makefile \
	   $(INSTALL_PATH)/examples/moduleGuide/scalarTransport/MTR
	cp $(wildcard scalarTransport/tracers/*.vars \
	              scalarTransport/tracers/*.vog \
	              scalarTransport/tracers/*.loci \
	              scalarTransport/tracers/*.lh \
	              scalarTransport/tracers/*.h) Makefile \
	   $(INSTALL_PATH)/examples/moduleGuide/scalarTransport/tracers

moduleGuide: FRC
	INSTALL_PATH="$(INSTALL_PATH)" $(MAKE) -C moduleGuide install

%/.: FRC
	mkdir -p  $(INSTALL_PATH)/examples/moduleGuide/$(firstword $(subst /, , $*))
	cp $(wildcard $(firstword $(subst /, , $*))/*.vog \
	   $(firstword $(subst /, , $*))/*.vars \
	   $(firstword $(subst /, , $*))/*.loci \
	   $(firstword $(subst /, , $*))/Make* \
           ) $(INSTALL_PATH)/examples/moduleGuide/$(firstword $(subst /, , $*))


FRC:
