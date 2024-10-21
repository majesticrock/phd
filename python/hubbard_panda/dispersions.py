import numpy as np
import matplotlib.pyplot as plt
import __path_appender as __ap
__ap.append()

from get_data import load_panda, hubbard_params
pd_data = load_panda("hubbard/square", "test", "dispersions.json.gz", **hubbard_params(0.0, -2.5, 0.0))

import dispersions_2D as d2d
import plot_settings as ps

# Initialize resolvents
resolvents = d2d.Dispersions2D(pd_data, messages=False)

# Define frequency range
w_lin = np.linspace(0., pd_data["continuum_boundaries"][1] + 0.3, 5000, dtype=complex)
w_lin += 1e-4j

# Set up the grid for the heatmap
N_index = len(pd_data["resolvents.amplitude_SC_a"])
index_range = np.arange(0, N_index + 1)  # Change this range as needed
omega_range = np.real(w_lin)

# Initialize a matrix to hold spectral densities
spectral_matrix_phase = np.zeros((len(w_lin), len(index_range)))
spectral_matrix_amplitude = np.zeros((len(w_lin), len(index_range)))

# Loop over indices and fill the matrix
for i, idx in enumerate(index_range - 1):
    try:
        spectral_matrix_phase[:, i]     = np.clip(np.real(resolvents.spectral_density(w_lin, "phase_SC_a", index=idx%N_index)), 0, 1)
        spectral_matrix_amplitude[:, i] = np.clip(np.real(resolvents.spectral_density(w_lin, "amplitude_SC_a", index=idx%N_index)), 0, 1)
    except TypeError:
        # we dont care
        continue

# Plot the phase part as a heatmap
fig, ax = plt.subplots(1, 2, figsize=(12, 6))

# Plot heatmap for phase
c_phase = ax[0].contourf(index_range, omega_range, spectral_matrix_phase, levels=100, cmap="viridis")
ax[0].set_ylabel(r"$\omega [t]$")
ax[0].set_xlabel(r"$\vec k$")
ax[0].set_title("Phase SC")
cbar_phase = fig.colorbar(c_phase, ax=ax[0])
cbar_phase.set_label(r"$A_\mathrm{Phase} (\omega) [t^{-1}]$")

# Plot heatmap for amplitude (Higgs) part
c_amplitude = ax[1].contourf(index_range, omega_range, spectral_matrix_amplitude, levels=100, cmap="inferno")
ax[1].set_ylabel(r"$\omega [t]$")
ax[1].set_xlabel(r"$\vec k$")
ax[1].set_title("Amplitude SC (Higgs)")
cbar_amp = fig.colorbar(c_amplitude, ax=ax[1])
cbar_amp.set_label(r"$A_\mathrm{Higgs} (\omega) [t^{-1}]$")

xticks_positions = [0, N_index // 3, 2 * N_index // 3, N_index]
xticks_labels = [r"$\Gamma$", "X", "M", r"$\Gamma$"]

# Set custom ticks and labels for both subplots
for axis in ax:
    axis.set_xticks(xticks_positions)
    axis.set_xticklabels(xticks_labels)

fig.tight_layout()
plt.show()
