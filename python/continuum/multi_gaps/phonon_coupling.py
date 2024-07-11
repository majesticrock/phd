import matplotlib.pyplot as plt
import numpy as np

import matplotlib.cm as cm
from matplotlib.colors import ListedColormap

import __path_appender as __ap
__ap.append()
from get_data import *

MEV_FACTOR = 1e3

TYPE = "g"
LABEL= r"$g$"
X_BOUNDS = [-0.1, 0.1]

main_df = load_all("continuum/exact_2000", "gap.json.gz").query(
            'omega_D == 10 & coulomb_scaling == 1 & E_F == 9.3')
main_df.sort_values(TYPE, inplace=True)
main_df.reset_index(inplace=True)

colors = cm.viridis(np.linspace(0, 1, len(main_df.index)))
my_cmap = ListedColormap(colors[:,:-1])
fig, ax = plt.subplots(1, 1)
sm = plt.cm.ScalarMappable(cmap=my_cmap, 
                norm=plt.Normalize(vmin=main_df[TYPE].min(), vmax=main_df[TYPE].max()))

for index, pd_row in main_df.iterrows():
    pd_row["data"]["ks"] -= pd_row["k_F"]
    plot_range = pd_row["data"].query(f'ks > {1.2 * X_BOUNDS[0]} & ks < {1.2 * X_BOUNDS[1]}')
    plot_data = (plot_range["Delta_Phonon"] + plot_range["Delta_Coulomb"]).to_numpy()
    if plot_data[int(len(plot_data) / 2)] < 0:
        plot_data *= -1
    ax.plot(plot_range["ks"], plot_data, ls="-", color=colors[index])

cbar = fig.colorbar(sm, ax=ax, label=f'{LABEL} [eV]')
ax.set_xlim(*X_BOUNDS)
ax.set_xlabel(r"$k - k_\mathrm{F} [\sqrt{\mathrm{eV}}]$")
ax.set_ylabel(r"$\Delta [\mathrm{meV}]$")

fig.tight_layout()

plt.show()