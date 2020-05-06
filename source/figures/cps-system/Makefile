FIGPATH  := ../../../figures
FIGNAME := $(notdir $(shell pwd))
KEYEXPORT := ../../../tools/keynote-export.applescript

.PHONY: clean cleanall

all: $(FIGNAME).pdf
	@cp $^ $(FIGPATH)

$(FIGNAME).pdf: $(FIGNAME).key
	@osascript $(KEYEXPORT) $(shell pwd)/$^ $(shell pwd)/$@ && pdfcrop $@ $@
	
clean:
 
cleanall:
