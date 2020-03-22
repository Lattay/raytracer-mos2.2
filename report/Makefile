LATEX=lualatex
FLAGS=

SRC=$(wildcard *.tex)
BASENAME=rapport
PDF=$(BASENAME).pdf
PDFVIEWER=zathura

.PHONY: clean main

main: $(PDF)

$(PDF): $(SRC) | build
	$(LATEX) -output-directory=build -jobname=$(BASENAME) $^
	mv build/$(PDF) $(PDF)

show:
	($(PDFVIEWER) $(PDF) &)

build:
	mkdir build

clean:
	rm -f *.log
	rm -rf build
