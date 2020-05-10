FIGPATH := ../../../figures
FIGNAME := $(notdir $(shell pwd))
LATEX   := latexmk -pdf -silent
PYTHON  := python3

.PHONY: clean cleanall

all: $(FIGNAME).pdf
	@cp $^ $(FIGPATH)

$(FIGNAME).pdf: $(FIGNAME).tex
	@$(LATEX) $(FIGNAME)

clean:
	@latexmk -c -pdf -quiet $(FIGNAME)

cleanall:
	@latexmk -C -pdf -quiet $(FIGNAME)
