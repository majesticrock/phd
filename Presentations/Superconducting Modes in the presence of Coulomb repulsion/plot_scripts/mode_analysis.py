import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import __path_appender as __ap
__ap.append()

## Testing setup
XTYPE= "g"

df = pd.read_pickle("modes/higgs_0.pkl").sort_values("g")
filtered_df = df[df['energies'].apply(len) > 0].reset_index()
fig, ax = plt.subplots()

for i in range(4):
    first_entries = filtered_df["energies"].apply(lambda x: x[i] if len(x) > i else None)
    valid_idx = first_entries.first_valid_index()  # Get the first index where first_entries is not None

    if valid_idx is not None:  # Ensure there's a valid entry to plot
        ax.plot(filtered_df[XTYPE], first_entries, color=f"C{i}")
        ax.scatter(filtered_df[XTYPE].iloc[valid_idx], first_entries.iloc[valid_idx], color=f"C{i}", marker="*")
        ax.axvline(filtered_df[XTYPE].iloc[valid_idx], ymax=first_entries.iloc[valid_idx], color=f"C{i}", ls="--")

ax.set_xlabel(r"$\Delta_\mathrm{max}$ $[\mathrm{meV}]$")
ax.set_ylabel(r"$\omega_0$ $[\mathrm{meV}]$")
ax.set_xlim(0, filtered_df[XTYPE].max())

fig.tight_layout()
fig.savefig(f"plots/{os.path.basename(__file__).split('.')[0]}_test.pdf")

N_MODES = 4
N_SETS = 3
scatter_data = []

MU_STAR = [0, 0.0526, 0.291]

fig, ax = plt.subplots()
for j in range(N_SETS):
    for spectral_type, marker in zip(["higgs", "phase"], ["*", "^"]):
        df = pd.read_pickle(f"modes/{spectral_type}_{j}.pkl").sort_values("g")
        filtered_df = df[df['energies'].apply(len) > 0].reset_index()

        for i in range(N_MODES):
            first_entries = filtered_df["energies"].apply(lambda x: x[i] if len(x) > i else None)
            valid_idx = first_entries.first_valid_index()
            if valid_idx is not None:
                x = filtered_df[XTYPE].iloc[valid_idx]
                if XTYPE == "g":
                    x -= MU_STAR[j]
                y = first_entries.iloc[valid_idx]
                scatter_data.append([x, y, j])
                
                ax.plot(x, y, color=f"C{j}", linestyle=None, marker=marker, markersize=10)

scatter_data = np.array(scatter_data).transpose()
for i in range(len(scatter_data[0])):
    print(scatter_data[0][i])
from ez_fit import ez_linear_fit
popt, pcov = ez_linear_fit(scatter_data[0], scatter_data[1], ax, x_bounds=[0, 1.1 * scatter_data[0].max()], zorder=-20, color="k", ls="--")

ax.text(0.02, 0.57, f"$a = ({popt[0]:1.3f} \\pm {np.sqrt(pcov[0][0]):1.3f})$", transform=ax.transAxes)
ax.text(0.02, 0.47, f"$b = ({popt[1]:1.3f} \\pm {np.sqrt(pcov[1][1]):1.3f})\\mathrm{{meV}}$", transform=ax.transAxes)

from matplotlib.lines import Line2D
# Create legend for marker styles
legend_markers = [Line2D([0], [0], marker='*', color='w', markerfacecolor='k', markersize=10, label=r'From $\mathcal{A}_\mathrm{Higgs}$'),
                  Line2D([0], [0], marker='^', color='w', markerfacecolor='k', markersize=10, label=r'From $\mathcal{A}_\mathrm{Phase}$'),
                  Line2D([0], [0], linestyle="--", color='k', label='Fit')]

# Create legend for colors
legend_colors = [Line2D([0], [0], color='C0', label='No Coulomb'),
                 Line2D([0], [0], color='C1', label='$\\lambda=1$'),
                 Line2D([0], [0], color='C2', label='$\\lambda=10^{-4}$')]

# Add the legends to the plot
old_leg = ax.legend(handles=legend_markers, loc='upper left')
ax.legend(handles=legend_colors, loc='lower right')
ax.add_artist(old_leg)

ax.set_xlabel(r"$\Delta_\mathrm{max}$ $[\mathrm{meV}]$")
ax.set_ylabel(r"$\omega_0$ $[\mathrm{meV}]$")

ax.set_xlim(0, 1.05 * scatter_data[0].max())
ax.set_ylim(0, 1.05 * scatter_data[1].max())

fig.tight_layout()
fig.savefig(f"plots/{os.path.basename(__file__).split('.')[0]}.pdf")
