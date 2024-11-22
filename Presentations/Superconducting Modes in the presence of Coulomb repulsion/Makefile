all: build/main.pdf

PLTS_SCR=gap_theta.py gap_coulomb.py resolvent_overview.py gap_behavior.py phase_peaks.py discretization.py \
	fit_gap_dependence.py energy_plot.py heatmap_resolvent.py max_true.py mode_analysis.py

#phase_limits.py peak_fits.py higgs_peak.py 

PLTS_PDF=$(addprefix plots/, $(subst .py,.pdf,$(PLTS_SCR)))

TeXOptions = -lualatex \
			 -interaction=nonstopmode \
			 -halt-on-error \
			 -output-directory=build

build/main.pdf: $(PLTS_PDF)

build/main.pdf: FORCE | build
	latexmk $(TeXOptions) main.tex 1> build/log || cat build/log

plots/%.pdf: plot_scripts/%.py matplotlibrc | plots
	TEXINPUTS=$$(pwd): python $<

modes/higgs_0.pkl: plots/heatmap_resolvent.pdf | modes

modes/all_data.pkl: plot_scripts/__all_data_pickler.py | modes
	python $<

plots/mode_analysis.pdf: modes/higgs_0.pkl
plots/heatmap_resolvent.pdf: modes/all_data.pkl
plots/max_true.pdf: plot_scripts/max_true.py modes/all_data.pkl
	TEXINPUTS=$$(pwd): python $<

FORCE:

modes:
	mkdir -p modes/

plots:
	mkdir -p plots/

build:
	mkdir -p build/

clean:
	rm -rf build/
	rm -rf plots/
	rm -rf modes/