MAIN      := thesis
DATADIR   := $(shell pwd)
DEFAULTS  := defaults
TEMPLATES := templates
SOURCE    := source
PARTS     := parts
BUILD     := build


PANDOC    := pandoc --data-dir $(DATADIR)
LATEX     := latexmk --outdir=$(BUILD) -pdf -use-make

BIB		  := $(HOME)/library.bib
CSL       := $(DATADIR)/aps.csl

MDS       := $(wildcard $(SOURCE)/*.md)
MDS2      := $(patsubst $(SOURCE)/%,$(PARTS)/%,$(MDS))


all: $(MAIN).pdf

$(MAIN).pdf: $(MAIN).tex
	@$(LATEX) $< 
	@cp $(BUILD)/$@ $@

$(MAIN).tex: $(MDS2) config.yaml $(TEMPLATES)/$(MAIN).latex $(DEFAULTS)/$(MAIN).yaml
	@$(PANDOC) \
			   --defaults $(MAIN) \
			   --output=$@ \
			   --bibliography=$(BIB)
			   

$(PARTS)/%.md: $(SOURCE)/%.md $(TEMPLATES)/%.markdown config.yaml | $(PARTS)
	@$(PANDOC) --from=markdown --to=gfm \
			   --template=$*.markdown \
			   --output=$@ \
			   --metadata-file=config.yaml \
			   --defaults empty-input \
			   --include-after $<

$(PARTS):
	@mkdir -p $@

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

