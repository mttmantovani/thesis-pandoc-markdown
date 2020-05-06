MAIN      := thesis
DATADIR   := $(shell pwd)
DEFAULTS  := defaults
TEMPLATES := templates
SOURCE    := source
PARTS     := parts
BUILD     := build


PANDOC    := pandoc --data-dir $(DATADIR)
LATEX     := latexmk --outdir=$(BUILD) -pdf -use-make -silent

BIB		  := $(HOME)/library.bib
CSL       := $(DATADIR)/aps.csl

MDS       := $(wildcard $(SOURCE)/*.md)
MDS2      := $(patsubst $(SOURCE)/%,$(PARTS)/%,$(MDS))

FIGSRC      := $(SOURCE)/figures
FIGS 		:= $(shell find $(FIGSRC) -mindepth 1 -maxdepth 1 -type d)
FIGURES		:= $(addsuffix .pdf, $(notdir $(FIGS)))
EXTS		:= py tex tikz tikzstyle key


all: $(MAIN).pdf

$(MAIN).pdf: $(MAIN).tex $(BIB) $(FIGURES)
	@$(LATEX) $< 
	@cp $(BUILD)/$@ $@

$(MAIN).tex: $(MDS2) config.yaml $(TEMPLATES)/$(MAIN).latex $(DEFAULTS)/$(MAIN).yaml
	@$(PANDOC) --defaults $(MAIN) --output=$@ --bibliography=$(BIB) \
	-M autoEqnLabels
			   
$(PARTS)/%.md: $(SOURCE)/%.md $(TEMPLATES)/%.markdown config.yaml | $(PARTS)
	@$(PANDOC) --defaults $(PARTS) --template=$*.markdown --output=$@ --include-after $<

$(PARTS):
	@mkdir -p $@

.SECONDEXPANSION:
$(FIGURES): %.pdf: $$(foreach ext,$$(EXTS),$$(wildcard $(FIGSRC)/%/*.$$(ext)))
	@cd $(FIGSRC)/$* && make && make cleanall

exportbib: $(BUILD)/$(MAIN).bcf
	@biber --output-format=bibtex --output-resolve --output-fieldcase=lower \
		   --output-directory=$(DATADIR) --output-file=references.bib -w -q \
		   $< 
clean:
	latexmk -c $(MAIN)
	rm -f $(PARTS)/*

cleanall:
	latexmk -C $(MAIN)
	rm -f *.pdf 
	rm -rf $(PARTS) $(BUILD)

