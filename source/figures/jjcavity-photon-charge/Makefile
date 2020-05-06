FIGPATH := ../../../figures
FIGNAME := $(notdir $(shell pwd))
LATEX   := latexmk -xelatex -silent
DATA    := $(shell ls *.csv | sort -f)
PYTHON  := python3

.PHONY: clean cleanall

all: $(FIGNAME).pdf
	@cp $^ $(FIGPATH)

$(FIGNAME).pdf: $(FIGNAME).tex $(FIGNAME).tikz $(FIGNAME)-plot.pdf
	@$(LATEX) $< && pdfcrop $@ $@

$(FIGNAME)-plot.pdf: $(FIGNAME)-plot.tex $(FIGNAME)-plot.pgf
	@$(LATEX) $< && pdfcrop $@ $@

$(FIGNAME)-plot.pgf: $(FIGNAME)-plot.py $(DATA)
	@$(PYTHON) $^

clean:
	@latexmk -c -pdf -quiet $(FIGNAME)

cleanall:
	@latexmk -C -pdf -quiet $(FIGNAME)
	@rm -f $(FIGNAME).pgf
