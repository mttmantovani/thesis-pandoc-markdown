FIGPATH  := ../../../figures
FIGNAME := $(notdir $(shell pwd))
TEXRULE := latexmk -xelatex -silent
DATA    := $(shell ls *.csv | sort -f)
PYTHON  := python3

.PHONY: clean cleanall

all: $(FIGNAME).pdf
	@cp $^ $(FIGPATH)

$(FIGNAME).pdf: $(FIGNAME).tex $(FIGNAME).pgf
	@$(TEXRULE) $< && pdfcrop $@ $@

$(FIGNAME).pgf: $(FIGNAME).py $(DATA)
	@$(PYTHON) $^

$(FIGNAME).pdf: $(FIGNAME).tex $(FIGNAME).py $(DATA)

clean:
	@latexmk -c -xelatex -quiet  $(FIGNAME)

cleanall:
	@latexmk -C -xelatex -quiet $(FIGNAME)
	@rm -f $(FIGNAME).pgf
