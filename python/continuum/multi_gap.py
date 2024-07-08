from turtle import pd
import matplotlib.pyplot as plt
import numpy as np

import matplotlib.cm as cm
from matplotlib.colors import ListedColormap

import __path_appender as __ap
__ap.append()
from create_zoom import create_zoom
from get_data import *

MEV_FACTOR = 1e3

TYPE = "omega_D"
LABLE= r"$\omega_\mathrm{D}$"

pd_data = load_all("continuum/disc_2000", "gap.json.gz").query('g == 10')
pd_data.sort_values(TYPE, inplace=True)
pd_data.reset_index(inplace=True)

colors = cm.gist_rainbow(np.linspace(0, 1, len(pd_data.index)))
my_cmap = ListedColormap(colors[:,:-1])
fig, ax = plt.subplots(1, 1)
sm = plt.cm.ScalarMappable(cmap=my_cmap, 
                norm=plt.Normalize(vmin=pd_data[TYPE].min(), vmax=pd_data[TYPE].max()))

for index, pd_row in pd_data.iterrows():
    pd_row["ks"] -= pd_row["k_F"]
    plot_data = pd_row["Delta_Phonon"] + pd_row["Delta_Coulomb"]
    ax.plot(pd_row["ks"], plot_data, ls="-", color=colors[index % 20])
    print(index)
    
cbar = fig.colorbar(sm, ax=ax, label=f'{LABLE} [meV]')
max_arg = pd_data[TYPE].idxmax()
axins = create_zoom(ax, 0.5, 0.3, 0.45, 0.6, xlim=(-0.15, 0.15), 
            ylim=(1.1 * np.min(plot_data), 1.05 * np.max(plot_data)))

ax.set_xlabel(r"$k - k_\mathrm{F} [\sqrt{\mathrm{eV}}]$")
ax.set_ylabel(r"$\Delta [\mathrm{meV}]$")

import matplotlib.ticker as ticker
ticks_y = ticker.FuncFormatter(lambda x, pos: '{0:g}'.format(x*MEV_FACTOR))
ax.yaxis.set_major_formatter(ticks_y)
axins.yaxis.set_major_formatter(ticks_y)

fig.tight_layout()

plt.show()