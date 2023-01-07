TEST_FILES=adiabatic.test Twall.test lid.test periodicSA.test periodicBSL.test periodicSST.test periodicDHRL.test periodicWilcox98.test periodicWilcox08.test periodicRKE.test wallLawAdiabatic.test wallLawTwall.test


include $(FLOWPSI_BASE)/test.conf

default: TestResults

TestResults: $(TEST_FILES) FRC
	cat $(TEST_FILES) > TestResults; rm -f $(TEST_FILES)

FRC:

%.test : %.vars
	FLOWPSI="$(FLOWPSI)" EXTRACT="$(EXTRACT)" NDIFF="$(NDIFF)" ./tests.sh $*.vars

clean:
	rm -fr TEST_* TestResults $(TEST_FILES)


