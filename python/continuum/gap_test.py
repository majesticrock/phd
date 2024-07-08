from turtle import pd
import matplotlib.pyplot as plt
import numpy as np

import __path_appender as __ap
__ap.append()
from create_zoom import create_zoom
from get_data import *

MEV_FACTOR = 1e3

fig, ax = plt.subplots()

pd_data = load_panda("continuum", "test", "gap.json.gz", 
                     **continuum_params(0.0, 1.0, 9.3, 10., 0.01)).iloc[0]
pd_data["ks"] -= pd_data["k_F"]

#0.0020500367846213573
ax.plot(pd_data["ks"], pd_data["Delta_Phonon"] + pd_data["Delta_Coulomb"], "k-", label=r"$\Delta_\mathrm{SC}$")
ax.plot(pd_data["ks"], pd_data["Delta_Phonon"], "--", label=r"$\Delta_\mathrm{Phonon}$")
ax.plot(pd_data["ks"], pd_data["Delta_Coulomb"], "--", label=r"$\Delta_\mathrm{Coulomb}$")
#ax.plot(pd_data["ks"], pd_data["Delta_cut"], ":", label=r"$\Delta_\mathrm{cut}$")
ax.plot(pd_data["ks"], pd_data["Delta_Fock"], "-", label=r"$\Delta_\mathrm{Fock}$")

axins = create_zoom(ax, 0.1, 0.3, 0.3, 0.6, xlim=(-0.02, 0.02), ylim=(1.1 * np.min(pd_data["Delta_Coulomb"]), 1.05 * np.max(pd_data["Delta_Phonon"])))
import matplotlib.ticker as ticker
ticks_y = ticker.FuncFormatter(lambda x, pos: '{0:g}'.format(x*MEV_FACTOR))
ax.yaxis.set_major_formatter(ticks_y)
axins.yaxis.set_major_formatter(ticks_y)

ax.set_xlabel(r"$k - k_\mathrm{F} [\sqrt{\mathrm{eV}}]$")
ax.set_ylabel(r"$\Delta [\mathrm{meV}]$")

ax.legend()
fig.tight_layout()

plt.show()