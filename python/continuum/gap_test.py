import matplotlib.pyplot as plt
import numpy as np

import __path_appender as __ap
__ap.append()
from create_zoom import create_zoom
from get_data import *

X_BOUNDS = [-0.1, 0.1]

fig, ax = plt.subplots()

main_df = load_panda("continuum", "exact_2000", "gap.json.gz", 
                     **continuum_params(0.0, 1., 9.3, 10., 10.))
pd_data = main_df["data"]
pd_data["ks"] -= main_df["k_F"]

ax.plot(pd_data["ks"], pd_data["Delta_Phonon"] + pd_data["Delta_Coulomb"], "k-", label=r"$\Delta_\mathrm{SC}$")
ax.plot(pd_data["ks"], pd_data["Delta_Phonon"], "--", label=r"$\Delta_\mathrm{Phonon}$")
ax.plot(pd_data["ks"], pd_data["Delta_Coulomb"], "--", label=r"$\Delta_\mathrm{Coulomb}$")
ax.plot(pd_data["ks"], pd_data["Delta_Fock"], "-", label=r"$\Delta_\mathrm{Fock}$")

axins = create_zoom(ax, 0.4, 0.3, 0.3, 0.6, xlim=(-0.02, 0.02), ylim=(1.1 * np.min(pd_data["Delta_Coulomb"]), 1.05 * np.max(pd_data["Delta_Phonon"])))

ax.set_xlabel(r"$k - k_\mathrm{F} [\sqrt{\mathrm{eV}}]$")
ax.set_ylabel(r"$\Delta [\mathrm{meV}]$")

ax.legend()
fig.tight_layout()

plt.show()