
guide.pdf: guide.tex 
	pdflatex guide
	bibtex guide
	pdflatex guide
	pdflatex guide

clean:
	rm -f *.aux *.ps *.bbl *.blg *.log *.toc *.dvi *.ps *~
distclean:
	rm -f *.aux *.ps *.bbl *.blg *.log *.toc *.dvi *.ps *.pdf *~
