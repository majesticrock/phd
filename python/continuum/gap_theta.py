import matplotlib.pyplot as plt
import numpy as np

import __path_appender as __ap
__ap.append()
from get_data import *

X_BOUNDS = [-0.05, 0.05]

fig, ax = plt.subplots()

main_df = load_panda("continuum", "exact_2000", "gap.json.gz", **continuum_params(0.0, 1., 9.3, 10., 10.))
pd_data = main_df["data"]
pd_data["ks"] -= main_df["k_F"]
plot_range = pd_data.query(f'ks > {1.2 * X_BOUNDS[0]} & ks < {1.2 * X_BOUNDS[1]}')
ax.plot(plot_range["ks"], plot_range["Delta_Phonon"] + plot_range["Delta_Coulomb"], "-", c="C0", label=r"$\Delta_\mathrm{SC}$")
ax.plot(plot_range["ks"], plot_range["Delta_Phonon"], "-.", c="C0", label=r"$\Delta_\mathrm{Phonon}$")
ax.plot(plot_range["ks"], plot_range["Delta_Coulomb"], "--", c="C0", label=r"$\Delta_\mathrm{Coulomb}$")

main_df = load_panda("continuum", "approx_2000", "gap.json.gz", **continuum_params(0.0, 1., 9.3, 10., 10.))
pd_data = main_df["data"]
pd_data["ks"] -= main_df["k_F"]
plot_range = pd_data.query(f'ks > {1.2 * X_BOUNDS[0]} & ks < {1.2 * X_BOUNDS[1]}')
ax.plot(plot_range["ks"], plot_range["Delta_Phonon"] + plot_range["Delta_Coulomb"], "-", c="C1", label=r"$\Delta_\mathrm{SC}$")
ax.plot(plot_range["ks"], plot_range["Delta_Phonon"], "-.", c="C1", label=r"$\Delta_\mathrm{Phonon}$")
ax.plot(plot_range["ks"], plot_range["Delta_Coulomb"], "--", c="C1", label=r"$\Delta_\mathrm{Coulomb}$")

ax.set_xlim(*X_BOUNDS)
ax.set_xlabel(r"$k - k_\mathrm{F} [\sqrt{\mathrm{eV}}]$")
ax.set_ylabel(r"$\Delta [\mathrm{meV}]$")

ax.legend()
fig.tight_layout()

plt.show()