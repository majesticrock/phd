import matplotlib.pyplot as plt
import __path_appender as __ap
__ap.append()
from get_data import *
from alpha_zip import *

X_BOUNDS = [1 - 0.003, 1 + 0.003]
Gs = [0.15, 0.5, 1.]

fig, axes = plt.subplots(nrows=len(Gs), sharex=True, figsize=(6.4, 6.4))

for ax, g, plot_label in alpha_label(axes, Gs):
    main_df = load_panda("continuum", "offset_10", "gap.json.gz", 
                        **continuum_params(N_k=20000, T=0, coulomb_scaling=0, screening=0, k_F=4.25, g=g, omega_D=10))
    pd_data = main_df["data"]
    pd_data["ks"] /= main_df["k_F"]
    plot_range = pd_data.query(f'ks > {0.8 * X_BOUNDS[0]} & ks < {1.2 * X_BOUNDS[1]}')
    plot_range.plot("ks", "Delta_Phonon", ax=ax, label=r"$\Delta_\mathrm{CUT}$", ls="-", c="C0", legend=False)

    main_df = load_panda("continuum", "theta_approx", "gap.json.gz", 
                        **continuum_params(N_k=8000, T=0, coulomb_scaling=0, screening=0, k_F=4.25, g=g, omega_D=10))
    pd_data = main_df["data"]
    pd_data["ks"] /= main_df["k_F"]
    plot_range = pd_data.query(f'Delta_Phonon > 0')
    plot_range.plot("ks", "Delta_Phonon", ax=ax, label=r"$\Delta_\mathrm{approx}$", ls="--", c="C1", legend=False)
    delta_edges = [plot_range["ks"].iloc[0], plot_range["ks"].iloc[-1]]
    ax.plot([X_BOUNDS[0], delta_edges[0], delta_edges[0]], [0, 0, plot_range["Delta_Phonon"].iloc[0]], color="C1", ls="--")
    ax.plot([X_BOUNDS[1], delta_edges[1], delta_edges[1]], [0, 0, plot_range["Delta_Phonon"].iloc[0]], color="C1", ls="--")
    ax.set_ylabel(r"$\Delta [\mathrm{meV}]$")
    
    ax.text(0.015, 0.84, f"{plot_label} $g={g}$", transform=ax.transAxes)
    
ax.set_xlim(*X_BOUNDS)
ax.set_xlabel(r"$k / k_\mathrm{F}$")
axes[0].legend()

fig.subplots_adjust(hspace=0)

import os
fig.savefig(f"plots/{os.path.basename(__file__).split('.')[0]}.pdf")
