import numpy as np
import matplotlib.pyplot as plt

import mrock.continued_fraction as cf
from mrock.get_data import *

# Create a DataLoader instance.
# By default, this will look for data in the package's default data directory
# or in the directory specified by the MROCK_DATA_DIR environment variable.
data_loader = DataLoader()

# System and sampling parameters for the lattice cut data.
SYSTEM = 'sc'
N = 16000
params = lattice_cut_params(
    N=N,
    g=1.,
    U=0.,
    E_F=0,
    omega_D=0.02
)

# Load resolvent data for the specified lattice/system and parameters.
main_df = data_loader.load_panda("lattice_cut", f"./{SYSTEM}", "resolvents.json.gz", **params)

# Build a ContinuedFraction evaluator and ignore noisy coefficients at both ends.
resolvents = cf.ContinuedFraction(main_df, ignore_first=200, ignore_last=280)
# Print the continuum edge estimate (gap) for reference.
print("Delta_true = ", resolvents.continuum_edges()[0])

# Prepare the figure and axis with labels.
fig, ax = plt.subplots()
ax.set_xlabel(r"$\omega / W$")
ax.set_ylabel(r"$\mathcal{A} (\omega) / W^{-1}$")

# Construct a positive frequency grid and add a small negative imaginary part for causality.
w_lin = np.linspace(0, .5 * main_df["continuum_boundaries"][1], 10000, dtype=complex)
w_lin += 1e-4j

# Evaluate phase and amplitude (Higgs) spectral densities on the grid.
A_phase = resolvents.spectral_density(w_lin, "phase_SC", with_terminator=True)
A_higgs = resolvents.spectral_density(w_lin, "amplitude_SC", with_terminator=True)

# Plot both spectra and set vertical limits.
ax.plot(w_lin.real, A_phase, label="Phase")
ax.plot(w_lin.real, A_higgs, label="Higgs")
ax.set_ylim(-0.05, 3.5)

# Shade the continuum region and finalize plot limits, legend, layout, then display.
resolvents.mark_continuum(ax)
ax.set_xlim(np.min(w_lin.real), np.max(w_lin.real))
ax.legend()
fig.tight_layout()
plt.show()