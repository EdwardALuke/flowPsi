TEST_FILES= steadyState.test timeAccurateKEC.test timeAccurateBDF2.test


include $(FLOWPSI_BASE)/test.conf

default: TestResults

TestResults: $(TEST_FILES) FRC
	cat $(TEST_FILES) > TestResults; rm -f $(TEST_FILES)

FRC:

%.test : %.vars
	FLOWPSI="$(FLOWPSI)" CHEMISTRY_DATABASE="$(FLOWPSI_BASE)" EXTRACT="$(EXTRACT)" NDIFF="$(NDIFF)" ./tests.sh $*.vars

clean:
	rm -fr TEST_* TestResults $(TEST_FILES)


