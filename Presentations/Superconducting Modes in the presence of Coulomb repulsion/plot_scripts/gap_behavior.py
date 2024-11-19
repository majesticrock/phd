import matplotlib.pyplot as plt
import numpy as np

import __path_appender as __ap
__ap.append()
from get_data import *
import os

fig, ax = plt.subplots()

all_data = load_all("continuum/offset_10/N_k=20000/T=0.0", "resolvents.json.gz").query("k_F == 4.25 & lambda_screening >= 1e-4")
plot_data = all_data.query("coulomb_scaling == 1 & omega_D == 10 & g in [0.3, 0.5, 0.7]").sort_values("lambda_screening")

linestyles = {0.3: '-', 0.5: '--', 0.7: ':'}
for g_value, group in plot_data.groupby("g"):
    group.plot(x="lambda_screening", y="Delta_max", ax=ax, label=f"g={g_value}", linestyle=linestyles.get(g_value, 'solid'), legend=True)

ax.set_xlabel(r"$\lambda$")
ax.set_ylabel(r"$\Delta_\mathrm{max}$ $[\mathrm{meV}]$")
ax.set_xscale('log')

ax.set_ylim(0, 4.5)

fig.tight_layout()
fig.savefig(f"plots/{os.path.basename(__file__).split('.')[0]}.pdf")
