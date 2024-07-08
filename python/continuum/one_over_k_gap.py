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
N = int(len(pd_data["ks"]) / 2) + 200
pd_data["ks"] = 1. / (pd_data["ks"][N:])**2

#ax.plot(pd_data["ks"], pd_data["Delta_Phonon"][N:] + pd_data["Delta_Coulomb"][N:], "k-", label=r"$\Delta_\mathrm{SC}$")
ax.plot(pd_data["ks"], pd_data["Delta_Coulomb"][N:], "x", label=r"$\Delta_\mathrm{Coulomb}$")

x_lin = np.linspace(30 * pd_data["ks"].min(), 0., 200)
ax.plot(x_lin, x_lin * pd_data["k_infinity_factor"], "-")

import matplotlib.ticker as ticker
ticks_y = ticker.FuncFormatter(lambda x, pos: '{0:g}'.format(x*MEV_FACTOR))
ax.yaxis.set_major_formatter(ticks_y)

ax.set_xlabel(r"$1 / k^2 [1 / \mathrm{eV}]$")
ax.set_ylabel(r"$\Delta [\mathrm{meV}]$")

ax.legend()
fig.tight_layout()

plt.show()