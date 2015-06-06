report.pdf: report.tex spconf.sty bibl_conf.bib $(shell find graphics)
	pdflatex report.tex > /dev/null
	bibtex report > /dev/null
	pdflatex report.tex > /dev/null
	pdflatex report.tex >  /dev/null
	
debug: report.tex spconf.sty $(shell find graphics)
	pdflatex report.tex
	bibtex report
	pdflatex report.tex
	pdflatex report.tex

.PHONY: all clean cleanall spell
	
all: report.pdf

clean:
	rm -rf *.bbl *.blg *.aux *.log *~ *.bak *.ps *.dvi *.log *.out *.tmp 

cleanall:
	rm -rf *.bbl *.blg *.aux *.log *~ *.bak *.ps *.dvi *.log *.pdf svnver.tex *.out *.tmp 

spell:
	ispell -f ispell.words -t *.tex
