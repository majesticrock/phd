import matplotlib.pyplot as plt
import numpy as np

import __path_appender as __ap
__ap.append()
from get_data import *
from scipy import integrate

X_BOUNDS = [-0.05, 0.05]

fig, ax = plt.subplots()

main_df = load_panda("continuum", "offset_10", "gap.json.gz", 
                    **continuum_params(N_k=8000, T=0.0, coulomb_scaling=1.0, screening=1e-4, k_F=4.25, g=5., omega_D=10.))
pd_data = main_df["data"]
pd_data["ks"] /= main_df["k_F"]
plot_range = pd_data.query(f'ks > {1.2 * X_BOUNDS[0]} & ks < {1.2 * X_BOUNDS[1]}')
ax.text(0.01, 0.8, fr'$\int \Delta_\mathrm{{exact}} \mathrm{{d}}k / k_\mathrm{{F}} = {integrate.trapz(plot_range["Delta_Phonon"], x=plot_range["ks"]) / main_df["k_F"]:.5f}$ meV', transform=ax.transAxes)
plot_range.plot("ks", "Delta_Phonon", ax=ax, label=r"$\Delta_\mathrm{exact}$", ls="-", c="C0")

main_df = load_panda("continuum", "offset_10", "gap.json.gz", 
                    **continuum_params(N_k=8000, T=0.0, coulomb_scaling=1.0, screening=1e-4, k_F=4.25, g=5., omega_D=10.))
pd_data = main_df["data"]
pd_data["ks"] /= main_df["k_F"]
plot_range = pd_data.query(f'ks > {1.2 * X_BOUNDS[0]} & ks < {1.2 * X_BOUNDS[1]}')
ax.text(0.01, 0.7, fr'$\int \Delta_\mathrm{{approx}} \mathrm{{d}}k / k_\mathrm{{F}} = {integrate.trapz(plot_range["Delta_Phonon"], x=plot_range["ks"]) / main_df["k_F"]:.5f}$ meV', transform=ax.transAxes)
plot_range.plot("ks", "Delta_Phonon", ax=ax, label=r"$\Delta_\mathrm{approx}$", ls="-", c="C1")

ax.set_xlim(*X_BOUNDS)
ax.set_xlabel(r"$k / k_\mathrm{F}$")
ax.set_ylabel(r"$\Delta [\mathrm{meV}]$")

ax.legend()
fig.tight_layout()

import os
plt.savefig(f"python/build/{os.path.basename(__file__).split('.')[0]}.svg")
plt.show()