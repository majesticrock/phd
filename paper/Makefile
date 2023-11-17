all: build/collective_excitations.pdf

TeXOptions = -lualatex \
			 -interaction=nonstopmode \
			 -halt-on-error \
			 -output-directory=build


build/collective_excitations.pdf: FORCE | build
	latexmk $(TeXOptions) collective_excitations.tex

FORCE:

build:
	mkdir -p build/

clean: