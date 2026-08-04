# Code for Collective Excitations in Correlated Quantum Materials

This repository contains the code used for the mean-field calculations, collective-mode analyses, and plotting workflows behind the thesis and publications listed below. The root of this project gives the basic structure and workflow, while the detailed instructions, and execution notes live in the corresponding subdirectories.

## Basic workflow

1. Build and install the shared core library in [PhdUtility](PhdUtility/README.md).
2. Generate the commutator-related prerequisites needed by the main simulations using [FermionCommute](cpp/FermionCommute/README.md).
3. Run the main calculations; see [cpp/ContinuumSystem](cpp/ContinuumSystem/README.md), [cpp/Hubbard](cpp/Hubbard/README.md), and [cpp/LatticeCUT](cpp/LatticeCUT/README.md).
4. Use the example scripts in [plot_examples](plot_examples/README.md) to visualize the generated data.


## What belongs where

- [PhdUtility](PhdUtility/README.md): the shared C++ library, build system, and installation logic.
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

- Collective excitations in competing phases in two and three dimensions, J. Althüser & G. S. Uhrig, Physical Review B 109, 205153 (2024), https://doi.org/10.1103/PhysRevB.109.205153
- Collective modes in superconductors including Coulomb repulsion, J. Althüser & G. S. Uhrig, SciPost Physics 19, 067 (2025), https://doi.org/10.21468/SciPostPhys.19.3.067
- Enhanced Superconductivity in Proximity to Peaks in Densities of States, J. Althüser, I. M. Eremin & G. S. Uhrig, (preprint) arXiv:2512.11451 (2025), https://doi.org/10.48550/arXiv.2512.11451
- Secondary Collective Excitations in Intermediate to Strong-Coupling Superconductors, J. Althüser & G. S. Uhrig, (preprint) arXiv:2605.20059 (2026), https://doi.org/10.48550/arXiv.2605.20059
