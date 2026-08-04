# Example plot scripts

These are short example plotting scripts that demonstrate common workflows with the code: loading precomputed data via `DataLoader`, evaluating continued-fraction resolvents, visualizing operator amplitudes, and performing a simple fit for the transition temperature. Each script is standalone and intended as an example; make sure mrock (and its dependencies) are installed and that the `DataLoader` can find the data (via the package data directory or the `MROCK_DATA_DIR` environment variable).

- continuum_delta.py
  Loads and plots the gap components Δ(k) from the continuum model as a function of k.

- continuum_bcs_comparison.py 
  Compares the gap Δ(k) from two continuum data sets (exact result and those that used the BCS approximation) as a function of k/k_F​ near the Fermi surface. Useful to visually assess the accuracy of the approximation and to mark the momentum region where the approximate gap is nonzero.

- continuum_resolvent.py  
  Loads the continuum system resolvent data and evaluates spectral densities for the superconducting phase and amplitude channels using the continued-fraction evaluator. Plots spectra in meV and shades the continuum region.

- latticecut_delta.py  
  Loads a lattice-cut dataset for a given system (bcc) and plots the superconducting gap Δ as a function of energy ε − μ. Useful to inspect the energy dependence of the gap on a single lattice cut.

- latticecut_operator_amplitudes.py  
  Loads full diagonalization results and uses FullDiagPurger to produce a 3-row figure showing amplitude (Higgs), occupation, and phase as a function of energy (shifted by the chemical potential). Demonstrates visualization of full-diagonalization output.

- latticecut_resolvent.py  
  Loads lattice-cut resolvents for a superconducting system and evaluates phase and amplitude spectral densities on a positive-frequency grid using the continued-fraction evaluator. Useful to compare Higgs vs. phase response and to locate continuum edges.

- latticecut_tc_single.py  
  Loads temperature-dependent maximum-gap data and fits Δ^2(T) linearly near the transition to extract the slope/intercept, from which the gap prefactor A and transition temperature Tc are obtained (with error propagation). Plots the data and the square-root of the linearized fit for a visual check.

- hubbard_resolvent.py
  Loads the Hubbard resolvent data and evaluates spectral densities for the superconducting phase and amplitude channels using the continued-fraction evaluator. Plots spectra in meV and shades the continuum region.

## How to run
- Execute each script with Python (e.g. `python continuum_delta.py`).
- Ensure required packages are installed (mrock, numpy, matplotlib, scipy, ...).
- If your data is not in the package default directory, set MROCK_DATA_DIR to point to your data folder.

These scripts are minimal examples; modify parameters (discretizations, interaction strengths, etc.) as needed.

Naturally, the data must exist to plot it.
Either take them from [here be TuDoData once I put all the data there] or generate them using the corresponding executable.