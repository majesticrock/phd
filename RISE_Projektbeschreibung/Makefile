all: build/RISE_Projektbeschreibung.pdf

TeXOptions = -pdflatex \
			 -interaction=nonstopmode \
			 -halt-on-error \
			 -output-directory=build
                                                                                
build/RISE_Projektbeschreibung.pdf: FORCE | build
		TEXTINPUTS=build: \
		BIBINPUTS=build: \
		max_print_line=1048576 \
	latexmk $(TeXOptions) RISE_Projektbeschreibung.tex 1> build/log || cat build/log

FORCE:

build:
	mkdir -p build/

clean:
	rm -rf build