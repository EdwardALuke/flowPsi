TEST_FILES=cylinder.test fixedMass.test inflow.test \
	isentropic.test nozzle.test superSonic.test wedge.test \
	shock_tube.test 


include $(FLOWPSI_BASE)/test.conf

default: TestResults

TestResults: $(TEST_FILES) FRC
	cat $(TEST_FILES) > TestResults; rm -f $(TEST_FILES)

FRC:

%.test : %.vars
	FLOWPSI="$(FLOWPSI)" EXTRACT="$(EXTRACT)" NDIFF="$(NDIFF)" ./tests.sh $*.vars

clean:
	rm -fr TEST_* TestResults $(TEST_FILES)


