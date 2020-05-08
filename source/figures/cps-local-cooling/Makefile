FIGPATH := ../../../figures
FIGNAME := $(notdir $(shell pwd))
LATEX   := latexmk -pdf -silent
DATA    := $(wildcard data/*)
PYTHON  := python3

.PHONY: clean cleanall

all: $(FIGNAME).pdf
	@cp $^ $(FIGPATH)

$(FIGNAME).pdf: $(FIGNAME).tex $(FIGNAME)-plot.pdf
	@$(LATEX) $(FIGNAME)

$(FIGNAME)-plot.pdf: $(FIGNAME)-plot.py $(DATA)
	@$(PYTHON) $^

clean:
	@latexmk -c -pdf -quiet $(FIGNAME)

cleanall:
	@latexmk -C -pdf -quiet $(FIGNAME)
	@rm -f $(FIGNAME)-plot.pdf
