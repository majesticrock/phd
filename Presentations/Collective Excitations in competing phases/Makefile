all: build/main.pdf

PLTS_SCR=phase_diagram.py resolvent_overview_SC_CDW.py resolvent_overview_AFM_CDW.py zero_peaks.py finite_size.py \
	sc_peak_in_cdw.py cdw_peak_in_sc.py resolvent_overview_V0.py afm_cdw_peaks_overview.py afm_cdw_peaks_details.py
	
#phase_limits.py peak_fits.py higgs_peak.py 

PLTS_PDF=$(addprefix plots/, $(subst .py,.pdf,$(PLTS_SCR)))

TeXOptions = -lualatex \
			 -interaction=nonstopmode \
			 -halt-on-error \
			 -output-directory=build

build/main.pdf: $(PLTS_PDF)

build/main.pdf: FORCE | build
	latexmk $(TeXOptions) main.tex 1> build/log || cat build/log

plots/%.pdf: %.py | plots
	TEXINPUTS=$$(pwd): python $<

FORCE:

plots:
	mkdir -p plots/

build:
	mkdir -p build/

clean:
	rm -rf build/
	rm -rf plots/