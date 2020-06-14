FIGPATH  := ../../../figures
FIGNAME := $(notdir $(shell pwd))
TEXRULE := latexmk -pdf -silent

.PHONY: clean cleanall

all: $(FIGNAME).pdf
	@cp $^ $(FIGPATH)

$(FIGNAME).pdf: $(FIGNAME).tex
	@$(TEXRULE) $<

clean:
	@latexmk -c -pdf -quiet $(FIGNAME)

cleanall:
	@latexmk -C -pdf -quiet $(FIGNAME)
