FIGPATH := ../../../figures
FIGNAME := $(notdir $(shell pwd))
LATEX   := latexmk -pdf -silent
DATA    := $(wildcard *.csv)
PYTHON  := python3

.PHONY: clean cleanall

all: $(FIGNAME).pdf
	@cp $^ $(FIGPATH)

$(FIGNAME).pdf: $(FIGNAME).py $(DATA)
	@$(PYTHON) $^

clean:

cleanall:
	@rm -f $(FIGNAME).pdf
