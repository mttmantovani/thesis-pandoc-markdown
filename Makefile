NAME      := thesis
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

# Build target directories if needed
$(info $(shell mkdir -p $(PARTS) $(BUILD)))

all: $(NAME).pdf

$(NAME).pdf: $(NAME).tex
	@$(LATEX) $< 
	@cp $(BUILD)/$@ $@

$(NAME).tex: $(MDS2) config.yaml $(TEMPLATES)/$(NAME).latex $(DEFAULTS)/$(NAME).yaml
	@$(PANDOC) \
			   --defaults thesis \
			   --output=$@ \
			   --bibliography=/Users/mattia/library.bib 
		# 	   config.yaml \
		# 	   parts/titlepage.md parts/preface.md parts/abstract.md \
		# 	   parts/zusammenfassung.md parts/toc.md parts/introduction.md \
		# 	   parts/theory.md parts/qdlaser.md parts/cps.md \
		# 	   parts/jjcavity.md \
		# 	   --bibliography=/Users/mattia/library.bib \
		# 	   --biblatex \
		# 	   --template=thesis.latex \
		# 	   --filter=pandoc-xnos
		# #	   parts/zusammenfassung.md parts/introduction.md \
		# #	   parts/theory.md parts/qdlaser.md parts/cps.md \
		# #	   parts/jjcavity.md \
			   

$(PARTS)/%.md: $(SOURCE)/%.md $(TEMPLATES)/%.markdown config.yaml
	@$(PANDOC) --from=markdown --to=gfm \
			   --template=$*.markdown \
			   --output=$@ \
			   --metadata-file=config.yaml \
			   --defaults empty-input \
			   --include-after $<



#thesis-print.pdf: thesis-print.tex 
#	@$(LATEX) $< && cp $(BUILD)/$@ $@
	
#thesis-print.tex: $(MDS) $(PARTS)/titlepage.tex $(DEFAULTS)/thesis-print.yaml
#	@$(PANDOC) --defaults thesis-print --variable print

#$(PARTS)/titlepage.tex: $(TEMPLATES)/titlepage.latex $(DEFAULTS)/titlepage.yaml metadata.yaml
#	@$(PANDOC) --defaults titlepage

clean:
	latexmk -c $(NAME)
	rm -f $(PARTS)/*

cleanall:
	latexmk -C $(NAME)
	rm -f *.pdf
	rm -rf $(PARTS) $(BUILD)

