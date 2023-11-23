all: build/main.pdf

TeXOptions = -lualatex \
			 -interaction=nonstopmode \
			 -halt-on-error \
			 -output-directory=build


build/collective_excitations.pdf: FORCE | build
	latexmk $(TeXOptions) collective_excitations.tex 1> build/log || cat build/log

build/main.pdf: FORCE | build
	latexmk $(TeXOptions) main.tex 1> build/log || cat build/log


FORCE:

build:
	mkdir -p build/

clean:
	rm -rf build/