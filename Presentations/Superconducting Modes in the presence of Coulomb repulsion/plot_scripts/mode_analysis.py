import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import __path_appender as __ap
__ap.append()

## Testing setup
XTYPE= "Delta_max"
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
        print(mode.x[0])



##### Mode gs | Deltas (meV):
# Higgs NC:
# 0.375  |  1.177819092419738
# 1.175  |  9.832135456722419
# 1.775  |  16.57786336321194
# 2.425  |  23.47804503837092
# Phase NC:
# 0.21   |  0.150247349509839
# 0.925  |  6.945049457183288
# 1.475  |  13.25736233398503
# 2.1    |  20.06692548862350
# 2.75   |  26.84033588256404
# 3.4    |  33.48488042464867
#
# Higgs lambda=1:
# 0.61   |  2.286570760619734
# 1.375  |  9.78876980708358
# 2.075  |  16.74739515262914
# 2.825  |  23.73280808430004
# Phase lambda=1
# 0.28   |  0.136015918311148
# 0.5    |  1.385479167671124
# 1.075  |  6.724192327275340
# 1.725  |  13.32830253769171
# 2.475  |  20.51726832933394
# 3.25   |  27.57445664499659
#
# Higgs lamda=1e-4:
# 1.1    |  2.881411527758234
# 2.0    |  9.368222013941743
# 2.95   |  16.25764521851743
# 3.9    |  22.58474283902128
# Phase lamda=1e-4:
# 1.625  |  6.566358359546582
# 2.475  |  12.88880493827046
# 3.425  |  19.47402483200409
# 

ax.set_xlabel(r"$\Delta_\mathrm{max}$ $[\mathrm{meV}]$")
ax.set_ylabel(r"$\omega_0$ $[\mathrm{meV}]$")
ax.set_xlim(0, filtered_df[XTYPE].max())

fig.tight_layout()
fig.savefig(f"plots/{os.path.basename(__file__).split('.')[0]}_test.pdf")


## Fitting setup
scatter_data = []
fig, ax = plt.subplots()
for j in range(N_SETS):
    for spectral_type, marker in zip(["higgs", "phase"], ["s", "o"]):
        df = pd.read_pickle(f"modes/{spectral_type}_{j}.pkl").sort_values("g")
        filtered_df = df[df['energies'].apply(len) > 0].reset_index()

        for i in range(N_MODES):
            if j < 2 and i == 0 and spectral_type=="phase": # standard phase mode is not relevant to us
                continue
            ith_entries = filtered_df["energies"].apply(lambda x: x[i] if len(x) > i else None)
            valid_idx = ith_entries.first_valid_index()
            if valid_idx is not None:
                y = filtered_df[XTYPE].iloc[valid_idx]
                x = i + int(j==2 and spectral_type=="phase") - 0.1 + 0.1 * j # ith_entries.iloc[valid_idx]
                scatter_data.append([ x, y ])
                # filtered_df[error] returns a list [lower, upper], but matplotlib requires an (2,n), i.e., (2,1) shaped object
                # Thus, we create a numpy array [[lower, upper]] and transpose it
                ax.errorbar(x, y, yerr=np.array([filtered_df[f"error_{XTYPE}"].iloc[valid_idx]]).transpose(), 
                            color=f"C{j}", linestyle=None, marker=marker, markersize=5)

scatter_data = np.array(scatter_data).transpose()
from matplotlib.lines import Line2D
# Create legend for marker styles
legend_markers = [Line2D([0], [0], marker="s", color='w', markerfacecolor='k', markersize=10, label=r'From $\mathcal{A}_\mathrm{Higgs}$'),
                  Line2D([0], [0], marker="o", color='w', markerfacecolor='k', markersize=10, label=r'From $\mathcal{A}_\mathrm{Phase}$')]

# Create legend for colors
legend_colors = [Line2D([0], [0], color='C0', label='No Coulomb'),
                 Line2D([0], [0], color='C1', label='$\\lambda=1$'),
                 Line2D([0], [0], color='C2', label='$\\lambda=10^{-4}$')]

# Add the legends to the plot
old_leg = ax.legend(handles=legend_markers, loc='upper left')
ax.legend(handles=legend_colors, loc='lower right')
ax.add_artist(old_leg)

ax.set_ylabel(r"$\Delta_\mathrm{max}$ $[\mathrm{meV}]$")
ax.set_xlabel(r"Mode number")

ax.set_xlim(-0.15, 4.05)
ax.set_ylim(0, 1.05 * scatter_data[1].max())

fig.tight_layout()
fig.savefig(f"plots/{os.path.basename(__file__).split('.')[0]}.pdf")
