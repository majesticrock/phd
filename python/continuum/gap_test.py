from turtle import pd
import matplotlib.pyplot as plt
import numpy as np
import gzip

import __path_appender as __ap
__ap.append()
from create_zoom import create_zoom
from get_data import load_panda, continuum_params

MEV_FACTOR = 1e3

fig, ax = plt.subplots()

pd_data = load_panda("continuum", "test", "gap.json.gz", 
                     **continuum_params(0.0, 1.0, 9.3, 12.71, 0.01)).iloc[0]
pd_data["ks"] -= pd_data["k_F"]

ax.plot(pd_data["ks"], pd_data["Delta_Phonon"] + pd_data["Delta_Coulomb"], "k-", label=r"$\Delta_\mathrm{SC}$")
ax.plot(pd_data["ks"], pd_data["Delta_Phonon"], "--", label=r"$\Delta_\mathrm{Phonon}$")
ax.plot(pd_data["ks"], pd_data["Delta_Coulomb"], "--", label=r"$\Delta_\mathrm{Coulomb}$")
ax.plot(pd_data["ks"], pd_data["Delta_Fock"], "-", label=r"$\Delta_\mathrm{Fock}$")

create_zoom(ax, 0.1, 0.3, 0.3, 0.6, xlim=(-0.02, 0.02), ylim=(1.1 * np.min(pd_data["Delta_Coulomb"]), 1.05 * np.max(pd_data["Delta_Phonon"])))

ax.set_xlabel(r"$k - k_\mathrm{F} [\sqrt{\mathrm{eV}}]$")
ax.set_ylabel(r"$\Delta [\mathrm{meV}]$")

import matplotlib.ticker as ticker
ticks_y = ticker.FuncFormatter(lambda x, pos: '{0:g}'.format(x*MEV_FACTOR))
ax.yaxis.set_major_formatter(ticks_y)

ax.legend()
fig.tight_layout()

plt.show()