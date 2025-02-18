
moduleGuide.pdf: moduleGuide.tex
	pdflatex moduleGuide
	bibtex moduleGuide
	pdflatex moduleGuide
	pdflatex moduleGuide

install: moduleGuide.pdf
	mkdir -p $(INSTALL_PATH)/docs
	cp moduleGuide.pdf $(INSTALL_PATH)/docs/moduleGuide.pdf

clean:
	rm -f *.aux *.ps *.bbl *.blg *.log *.toc *.dvi *.ps *~
distclean:
	rm -f *.aux *.ps *.bbl *.blg *.log *.toc *.dvi *.ps moduleGuide.pdf *~
