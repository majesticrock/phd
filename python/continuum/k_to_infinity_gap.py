from turtle import pd
import matplotlib.pyplot as plt
import numpy as np

import __path_appender as __ap
__ap.append()
from get_data import *

MEV_FACTOR = 1e3

fig, ax = plt.subplots()

main_df = load_panda("continuum", "limits", "gap.json.gz", **continuum_params(0.0, 1., 9.3, 10., 10.))
pd_data = main_df["data"]
pd_data["total"] = np.abs( pd_data["Delta_Coulomb"])
pd_data.plot("ks", "total", ax=ax, label=r"$\Delta_\mathrm{Coulomb}$")

x_lin = np.linspace(0.1 * pd_data["ks"].max(), 2. * pd_data["ks"].max(), 200)
y_inf = 1e3 * np.abs(main_df["k_infinity_factor"]) / x_lin**2
y_zero = 1e3 * np.abs(main_df["k_zero_factor"])

ax.plot(x_lin, y_inf, "--", c="C1", label=r"$k\to\infty$")
ax.axhline(y_zero, ls="--", c="C3", label=r"$k\to 0$")

ax.set_xlabel(r"$k [\sqrt{\mathrm{eV}}]$")
ax.set_ylabel(r"$|\Delta|[\mathrm{meV}]$")

ax.set_yscale("log")
ax.set_xscale("log")

ax.legend()
fig.tight_layout()

import os
plt.savefig(f"python/build/{os.path.basename(__file__).split('.')[0]}.svg")
plt.show()