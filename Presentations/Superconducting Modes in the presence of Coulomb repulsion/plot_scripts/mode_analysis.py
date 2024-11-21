import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import __path_appender as __ap
__ap.append()

## Testing setup
XTYPE= "true_gap"
N_MODES = 5
N_SETS = 3
MU_STAR = [0, 0.0526, 0.291]

class Mode:
    def __init__(self, first_x, first_energy):
        self.x = [first_x]
        self.energies = [first_energy]
    def append(self, new_x, new_energy):
        self.x.append(new_x)
        self.energies.append(new_energy)

class ModeCollector:
    def __init__(self, first_x, first_energies):
        self.modes = [ Mode(first_x, energy) for energy in first_energies if energy is not None]
    
    def append_new_energy(self, new_x, new_energies):
        EPS = 1
        for new_energy in new_energies:
            if new_energy is None:
                continue
            best_mode_idx = None
            best_energy_diff = 10000
            for j, mode in enumerate(self.modes):
                if np.abs(mode.energies[-1] - new_energy) < best_energy_diff:
                    best_energy_diff = np.abs(mode.energies[-1] - new_energy)
                    best_mode_idx = j
            if best_mode_idx is None or best_energy_diff > EPS:
                self.modes.append(Mode(new_x, new_energy))
            else:
                self.modes[best_mode_idx].append(new_x, new_energy)

fig, ax = plt.subplots()
for idx, spectral_type in enumerate(["higgs", "phase"]):
    df = pd.read_pickle(f"modes/{spectral_type}_2.pkl").sort_values("g")
    filtered_df = df[df['energies'].apply(len) > 0].reset_index()
    
    for i, df_row in filtered_df.iterrows():
        if i == 0:
            modes = ModeCollector(df_row[XTYPE], df_row["energies"])
        else:
            modes.append_new_energy(df_row[XTYPE], df_row["energies"])
    
    ax.plot(filtered_df[XTYPE], filtered_df["true_gap"], "k-", zorder=-100)
    for mode in modes.modes:
        ax.plot(mode.x, mode.energies, color=f"C{idx}")
        ax.scatter(mode.x[0], mode.energies[0], color=f"C{idx}", marker="*")

ax.set_xlabel(r"$\Delta_\mathrm{max}$ $[\mathrm{meV}]$")
ax.set_ylabel(r"$\omega_0$ $[\mathrm{meV}]$")
ax.set_xlim(0, filtered_df[XTYPE].max())

fig.tight_layout()
fig.savefig(f"plots/{os.path.basename(__file__).split('.')[0]}_test.pdf")


## Fitting setup
scatter_data = []
fig, ax = plt.subplots()
for j in range(N_SETS):
    for spectral_type, marker in zip(["higgs", "phase"], ["*", "^"]):
        df = pd.read_pickle(f"modes/{spectral_type}_{j}.pkl").sort_values("g")
        filtered_df = df[df['energies'].apply(len) > 0].reset_index()

        for i in range(N_MODES):
            if j < 2 and i == 0 and spectral_type=="phase": # standard phase mode is not relevant to us
                continue
            ith_entries = filtered_df["energies"].apply(lambda x: x[i] if len(x) > i else None)
            valid_idx = ith_entries.first_valid_index()
            if valid_idx is not None:
                x = filtered_df[XTYPE].iloc[valid_idx]
                y = i + int(j==2 and spectral_type=="phase")# ith_entries.iloc[valid_idx]
                scatter_data.append([x, y, j])
                
                ax.plot(x, y, color=f"C{j}", linestyle=None, marker=marker, markersize=10)

scatter_data = np.array(scatter_data).transpose()
from matplotlib.lines import Line2D
# Create legend for marker styles
legend_markers = [Line2D([0], [0], marker='*', color='w', markerfacecolor='k', markersize=10, label=r'From $\mathcal{A}_\mathrm{Higgs}$'),
                  Line2D([0], [0], marker='^', color='w', markerfacecolor='k', markersize=10, label=r'From $\mathcal{A}_\mathrm{Phase}$')]

# Create legend for colors
legend_colors = [Line2D([0], [0], color='C0', label='No Coulomb'),
                 Line2D([0], [0], color='C1', label='$\\lambda=1$'),
                 Line2D([0], [0], color='C2', label='$\\lambda=10^{-4}$')]

# Add the legends to the plot
old_leg = ax.legend(handles=legend_markers, loc='upper left')
ax.legend(handles=legend_colors, loc='lower right')
ax.add_artist(old_leg)

ax.set_xlabel(r"$\Delta_\mathrm{max}$ $[\mathrm{meV}]$")
ax.set_ylabel(r"Mode number")

ax.set_xlim(0, 1.05 * scatter_data[0].max())

fig.tight_layout()
fig.savefig(f"plots/{os.path.basename(__file__).split('.')[0]}.pdf")
