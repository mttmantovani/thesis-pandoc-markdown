FIGPATH  := ../../../figures
FIGNAME := $(notdir $(shell pwd))
TEXRULE := latexmk -xelatex -silent
PYTHON  := python3
DATA    := $(wildcard *.txt)

.PHONY: clean cleanall

all: $(FIGNAME).pdf
	@cp $^ $(FIGPATH)

$(FIGNAME).pdf: $(FIGNAME).py
	@$(PYTHON) $< $(DATA) && pdfcrop $@ $@

clean:

cleanall:
	@rm -f $(FIGNAME).pdf

