FIGPATH  := ../
FIGNAME := $(notdir $(shell pwd))
TEXRULE := latexmk -xelatex 

.PHONY: clean cleanall

all: $(FIGNAME).pdf
	@cp $^ $(FIGPATH)

$(FIGNAME).pdf: $(FIGNAME).tex
	@$(TEXRULE) $<

clean:
	@latexmk -c 

cleanall:
	@latexmk -C
