import numpy as np
import matplotlib.pyplot as plt

import mrock.FullDiagPurger as fdp
from mrock.get_data import *

# Create a DataLoader instance.
# By default, this will look for data in the package's default data directory
# or in the directory specified by the MROCK_DATA_DIR environment variable.
data_loader = DataLoader()

# Basic run parameters and directory for the system data.
SYSTEM = 'sc'
E_F = 0
OMEGA_D = 0.02
G = 3
DIR = f"./{SYSTEM}"
N = 16000

# Prepare a 3-row stacked figure for Higgs (amplitude), Occupation, and Phase.
fig, axes = plt.subplots(nrows=3, sharex=True)
fig.subplots_adjust(hspace=0)
axes[0].set_ylabel("Higgs")
axes[1].set_ylabel("Occupation")
axes[2].set_ylabel("Phase")
axes[-1].set_xlabel(r"$\varepsilon - \mu$")
axes[0].set_xlim(-0.25, 0.25)

# Build lattice cut parameters (sampling and physical params) and load full diagonalization data.
params = lattice_cut_params(
    N=N,
    g=G,
    U=0,
    E_F=E_F,
    omega_D=OMEGA_D
)
main_df = data_loader.load_panda(
    "lattice_cut",
    DIR,
    "full_diagonalization.json.gz",
    print_date=False,
    **params
)

# Create a FullDiagPurger to process/visualize the diagonalization results.
# The energy axis is shifted by the chemical potential.
purger = fdp.FullDiagPurger(main_df, np.linspace(-1, 1, N) - main_df["chemical_potential"])

# Plot amplitude (Higgs) and occupation on the first two axes, and phase on the third.
purger.plot_amplitude(axes[:2], combined_norm=True)
purger.plot_phase(axes[2], label="Result")

# Add zero lines for reference and display.
for ax in axes:
    ax.axhline(0, c="k", ls=":")

plt.show()