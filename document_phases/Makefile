all: build/thesis.pdf

PLTS_SCR=

PLTS_PDF=$(addprefix build/, $(subst .py,.pdf,$(PLTS_SCR)))

TeXOptions = -lualatex \
			 -interaction=nonstopmode \
			 -halt-on-error \
			 -output-directory=build
                                                                                
build/thesis.pdf: $(PLTS_PDF)

build/thesis.pdf: FORCE | build
	latexmk $(TeXOptions) thesis.tex 1> build/log || cat build/log
	
build/%.pdf: plot_scripts/plot_%.py
	TEXINPUTS=$$(pwd): python $<

FORCE:

build:
	mkdir -p build/

clean:
	rm -rf build
