import numpy as np
import matplotlib.pyplot as plt

import mrock.continued_fraction as cf
from mrock.get_data import *

# Create a DataLoader instance.
# By default, this will look for data in the package's default data directory
# or in the directory specified by the MROCK_DATA_DIR environment variable.
data_loader = DataLoader()
pd_data = data_loader.load_panda(
    "continuum",
    "offset_25",
    "resolvents.json.gz",
    **continuum_params(
        N_k=30000,
        T=0,
        coulomb_scaling=1,
        screening=1e-4,
        k_F=4.25,
        g=3.65,
        omega_D=10
    )
)

# Build a ContinuedFraction evaluator and ignore noisy coefficients at the edges.
resolvents = cf.ContinuedFraction(pd_data, ignore_first=80, ignore_last=90)
# Print an estimate of the true gap from the continuum edge (converted to meV).
print("Delta_true = ", 0.5e3 * resolvents.continuum_edges()[0])

# Prepare the figure and axis with axis labels and display limits.
fig, ax = plt.subplots()
ax.set_ylim(-0.01, 0.1)
ax.set_xlabel(r"$\omega [\mathrm{meV}]$")
ax.set_ylabel(r"$\mathcal{A} (\omega) [\mathrm{eV}^{-1}]$")

# Construct a real frequency grid (add small negative imaginary part for causality).
w_lin = np.linspace(
    -0.005 * pd_data["continuum_boundaries"][1],
    1.5 * pd_data["continuum_boundaries"][0],
    5000,
    dtype=complex
)
w_lin += 1e-6j

# Plot the phase and amplitude (Higgs) spectral densities.
# Multiply frequencies by 1e3 to convert to meV for the x-axis.
ax.plot(
    1e3 * w_lin.real,
    resolvents.spectral_density(w_lin, "phase_SC", with_terminator=True),
    label="Phase"
)
ax.plot(
    1e3 * w_lin.real,
    resolvents.spectral_density(w_lin, "amplitude_SC", with_terminator=True),
    label="Higgs"
)

# Shade the continuum region (scale factor 1e3 for meV units).
resolvents.mark_continuum(ax, 1e3)

# Finalize plot limits, legend and layout, then display.
ax.set_xlim(1e3 * np.min(w_lin.real), 1e3 * np.max(w_lin.real))
ax.legend()
fig.tight_layout()
plt.show()