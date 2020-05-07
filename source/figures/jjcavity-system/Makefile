FIGPATH := ../../../figures
FIGNAME := $(notdir $(shell pwd))
LATEX := latexmk -pdf -silent

.PHONY: clean cleanall

all: $(FIGNAME).pdf
	@cp $< $(FIGPATH)

$(FIGNAME).pdf: $(FIGNAME).tex $(FIGNAME).tikz
	@$(LATEX) $<

clean:
	@latexmk -c -pdf -quiet $(FIGNAME)

cleanall:
	@latexmk -C -pdf -quiet $(FIGNAME)
