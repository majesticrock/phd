all: build/main.pdf

PLTS_SCR=phase_diagram.py

PLTS_PDF=$(addprefix plots/, $(subst .py,.pdf,$(PLTS_SCR)))

TeXOptions = -lualatex \
			 -interaction=nonstopmode \
			 -halt-on-error \
			 -output-directory=build

build/main.pdf: $(PLTS_PDF)

build/main.pdf: FORCE | build
	latexmk $(TeXOptions) main.tex 1> build/log || cat build/log

plots/%.pdf: %.py
	TEXINPUTS=$$(pwd): python $<

FORCE:

build:
	mkdir -p build/

clean:
	rm -rf build/