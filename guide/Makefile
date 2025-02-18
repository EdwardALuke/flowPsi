
flowPsiGuide.pdf: guide.tex  axi_source_terms.tex generalized_grid.tex govern.tex grid_terms.tex inviscidFlux.tex jacobian.tex scalar_transport.tex space_integration.tex time_integration.tex turbulence.tex sst.tex wilcox98.tex spalartallmaras.tex baseline.tex compressibility_corrections.tex dhrl.tex
	pdflatex guide
	bibtex guide
	pdflatex guide
	pdflatex guide
	mv guide.pdf flowPsiGuide.pdf

install: flowPsiGuide.pdf
	mkdir -p $(INSTALL_PATH)/docs
	cp flowPsiGuide.pdf $(INSTALL_PATH)/docs/flowPsiGuide.pdf

clean:
	rm -f *.aux *.ps *.bbl *.blg *.log *.toc *.dvi *.ps *~
distclean:
	rm -f *.aux *.ps *.bbl *.blg *.log *.toc *.dvi *.ps flowPsiGuide.pdf *~
