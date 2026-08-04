# Code for Collective Excitations in Correlated Quantum Materials

Computing spectral functions is a demanding task. Here, the iterated equations of motion scheme is used to address the issue.
Employing the method involves calculating various commutators and the reduction of high-order expectation values to computable bilinear expressions using Wick's theorem.
This daunting task is handled by the [symbolic_operators](PhdUtility/symbolic_operators/README.md) sublibrary located in `PhdUtility`.

The algorithms derived in [Ref. 1](https://doi.org/10.1103/PhysRevB.109.205153) and [Ref. 4](https://doi.org/10.48550/arXiv.2605.20059) are implemented in the [iEoM](PhdUtility/iEoM/README.md) sublibrary also located in `PhdUtility`.

The sublibrary [utility](PhdUtility/utility/README.md) (also located in `PhdUtility`) contains various functionality that is reused throughout the projects.

The `cpp` directory contains the three applications that employ the aforementioned libraries to evaluate the spectra of collective excitations in different systems.

The `plot_examples` directory contains a few example scripts for plotting the results using matplotlib in python.

The computationally demanding parts of the code are written in C++.

The root of this project gives the basic structure and workflow, while the detailed instructions, and execution notes live in the corresponding subdirectories.

## Basic workflow

1. Build and install the shared core library in [PhdUtility](PhdUtility/README.md).
2. Generate the commutator-related prerequisites needed by the main simulations using [FermionCommute](cpp/FermionCommute/README.md).
3. Run the main calculations; see [cpp/ContinuumSystem](cpp/ContinuumSystem/README.md), [cpp/Hubbard](cpp/Hubbard/README.md), and [cpp/LatticeCUT](cpp/LatticeCUT/README.md).
4. Use the example scripts in [plot_examples](plot_examples/README.md) to visualize the generated data.


## What belongs where

- [PhdUtility](PhdUtility/README.md): the shared C++ library, build system, and installation logic.
- [cpp/FermionCommute](cpp/FermionCommute/README.md): Application using the [symbolic_operators](PhdUtility/symbolic_operators/README.md) sublibrary to evaluate commutators and expectation values.
- [cpp/Hubbard](cpp/Hubbard/README.md): mean-filed and collective-mode calculations on the half-filled extended Hubbard model.
- [cpp/ContinuumSystem](cpp/ContinuumSystem/README.md): using an effective interaction in a continuum-model.
- [cpp/LatticeCUT](cpp/LatticeCUT/README.md): using an effective interaction on lattices, including full-diagonalization workflows.
- [plot_examples](plot_examples/README.md): Python scripts that demonstrate how to load and plot the generated data.
- [data](data/): the main location for generated data used by the analysis and plotting scripts. Can of course be a soft link.

## Typical prerequisites

A typical run needs:

- a C++20-capable compiler,
- CMake and a standard build toolchain,
- Eigen, Boost, OpenMP, and nlohmann/json,
- Python with matplotlib, numpy, and scipy for the plotting examples.

The exact requirements and parameter conventions for each project are documented in the relevant subdirectory README.

## Tests

The entire project can be tested using `test_everything.sh`.

## Related publications

This superproject contains all code relevant to my doctoral thesis (not yet published) and the following publications:

1. Collective excitations in competing phases in two and three dimensions, J. Althüser & G. S. Uhrig, Physical Review B 109, 205153 (2024), https://doi.org/10.1103/PhysRevB.109.205153
1. Collective modes in superconductors including Coulomb repulsion, J. Althüser & G. S. Uhrig, SciPost Physics 19, 067 (2025), https://doi.org/10.21468/SciPostPhys.19.3.067
1. Enhanced Superconductivity in Proximity to Peaks in Densities of States, J. Althüser, I. M. Eremin & G. S. Uhrig, (preprint) arXiv:2512.11451 (2025), https://doi.org/10.48550/arXiv.2512.11451
1. Secondary Collective Excitations in Intermediate to Strong-Coupling Superconductors, J. Althüser & G. S. Uhrig, (preprint) arXiv:2605.20059 (2026), https://doi.org/10.48550/arXiv.2605.20059
