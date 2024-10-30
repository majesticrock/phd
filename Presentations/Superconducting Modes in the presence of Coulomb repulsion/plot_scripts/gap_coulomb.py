import matplotlib.pyplot as plt
import numpy as np

import __path_appender as __ap
__ap.append()
from get_data import *
from create_zoom import *
import os
from alpha_zip import alpha_label

X_BOUNDS = [-0.1, 0.1]

fig, axes = plt.subplots(nrows=2, sharex=True, figsize=(6.4, 6.4))
for ax, screening, label in alpha_label(axes, [1e-4, 1]):
    main_df = load_panda("continuum", "offset_20", "gap.json.gz", 
                        **continuum_params(N_k=20000, T=0, coulomb_scaling=1, screening=screening, k_F=4.25, g=0.8, omega_D=10))
    pd_data = main_df["data"]
    pd_data["ks"] /= main_df["k_F"]

    if pd_data["Delta_Coulomb"][0] > 0:
        pd_data["Delta_Phonon"] *= -1
        pd_data["Delta_Coulomb"] *= -1
    ax.plot(pd_data["ks"], pd_data["Delta_Phonon"] + pd_data["Delta_Coulomb"], "k-", label=r"$\Delta$")
    pd_data.plot(x="ks", y=["Delta_Phonon", "Delta_Coulomb", "Delta_Fock"], ax=ax, style=['--', '--', ':'], 
                label=[r"$\Delta_\mathrm{Ph}$", r"$\Delta_\mathrm{C}$", r"$\epsilon_\mathrm{C}$"], legend=False)

    axins = create_zoom(ax, 0.15, 0.35, 0.27, 0.55, xlim=(1-0.005, 1.005), 
                        ylim=(1.7 * np.min(pd_data["Delta_Phonon"] + pd_data["Delta_Coulomb"]), 1.1 * np.max(pd_data["Delta_Phonon"])), mark_inset=False)
    ax.set_ylabel(r"$\Delta [\mathrm{meV}]$")
    ax.text(0.015, 0.88, label, transform=ax.transAxes)

axes[0].legend()
axes[1].set_xlabel(r"$k / k_\mathrm{F}$")

fig.subplots_adjust(hspace=0)
fig.savefig(f"plots/{os.path.basename(__file__).split('.')[0]}.pdf")