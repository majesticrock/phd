all: build/main.pdf

PLTS_SCR=gap_theta.py gap_coulomb.py resolvent_overview.py heatmap_resolvent.py gap_behavior.py phase_peaks.py discretization.py \
	fit_gap_dependence.py

#phase_limits.py peak_fits.py higgs_peak.py 

PLTS_PDF=$(addprefix plots/, $(subst .py,.pdf,$(PLTS_SCR)))

TeXOptions = -lualatex \
			 -interaction=nonstopmode \
			 -halt-on-error \
			 -output-directory=build

build/main.pdf: $(PLTS_PDF) plots/mode_analysis.pdf

build/main.pdf: FORCE | build
	latexmk $(TeXOptions) main.tex 1> build/log || cat build/log

plots/%.pdf: plot_scripts/%.py matplotlibrc | plots
	TEXINPUTS=$$(pwd): python $<

modes/higgs.pkl: plots/heatmap_resolvent.pdf

plots/mode_analysis.pdf: plot_scripts/mode_analysis.py modes/higgs_0.pkl
	TEXINPUTS=$$(pwd): python $<

FORCE:

plots:
	mkdir -p plots/

build:
	mkdir -p build/

clean:
	rm -rf build/
	rm -rf plots/