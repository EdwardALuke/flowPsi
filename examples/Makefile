
EXAMPLE_CASES=$(wildcard */*.vars)

.phony: FRC setup copy

setup: FRC
	mkdir -p $(INSTALL_PATH)/examples

copy :$(EXAMPLE_CASES) setup

install: copy moduleExamplesInstall


%.vars: FRC
	mkdir -p  $(INSTALL_PATH)/examples/$(firstword $(subst /, , $*))
	cp $(wildcard $(firstword $(subst /, , $*))/*.vog \
	                $(firstword $(subst /, , $*))/*.vars \
	                $(firstword $(subst /, , $*))/*.dat \
	                $(firstword $(subst /, , $*))/*.bc \
	                $(firstword $(subst /, , $*))/*.grd \
	                $(firstword $(subst /, , $*))/*.cc \
	) $(INSTALL_PATH)/examples/$(firstword $(subst /, , $*))

moduleExamplesInstall: setup FRC
	INSTALL_PATH=$(INSTALL_PATH) $(MAKE) -C moduleExamples install

FRC:
