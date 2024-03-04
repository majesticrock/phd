all: poster.pdf

.phony: clean new

clean:
	rm poster.pdf

new:
	pdflatex poster.tex

poster.pdf: poster.tex
	pdflatex poster.tex