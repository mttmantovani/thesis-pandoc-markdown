FIGPATH  := ../../../figures
FIGNAME := $(notdir $(shell pwd))
TEXRULE := latexmk -xelatex -silent
PYTHON  := python3

.PHONY: clean cleanall

all: $(FIGNAME).pdf
	@cp $^ $(FIGPATH)

$(FIGNAME).pdf: $(FIGNAME).py
	@$(PYTHON) $< && pdfcrop $@ $@

clean:

cleanall:
	@rm -f $(FIGNAME).pdf
