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

#$(FIGNAME).pdf: $(FIGNAME).tex $(FIGNAME).pgf
#	@$(TEXRULE) $< && pdfcrop $@ $@
#
#$(FIGNAME).pgf: $(FIGNAME).py $(DATA)
#	@$(PYTHON) $^
#
#$(FIGNAME).pdf: $(FIGNAME).tex $(FIGNAME).py $(DATA)
#
#clean:
#	@latexmk -c $(FIGNAME)
#
#cleanall:
#	@latexmk -C $(FIGNAME)
#	@rm -f $(FIGNAME).pgf
