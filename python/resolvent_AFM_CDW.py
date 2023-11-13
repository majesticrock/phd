import numpy as np
import matplotlib.pyplot as plt
import lib.continued_fraction as cf
from lib.iterate_containers import iterate_containers
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, mark_inset
# Calculates the resolvent in w^2

prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']

Ts = np.array([0.])
Us = np.array([3.8])
Vs = np.array([1.])

use_XP = True
createZoom = False

folder = "data/modes/square/dos_900/"
name_suffix = "CDW"
element_names = ["a", "a+b", "a+ib"]
fig, ax = plt.subplots()

#ax.set_xscale("log")
ax.set_yscale("log")
#ax.set_ylim(0, 1)

plot_lower_lim = 2
plot_upper_lim = 10

if createZoom:
    axins = inset_axes(ax, width='50%', height='30%', loc='center right')
    # Plot the zoomed-in region
    axins.set_xlim(2.75, 2.85)
    axins.set_ylim(0, 1)
    plt.yticks(visible=False)
    plt.xticks(visible=False)
    # Mark the area in the main plot with a rectangle and a connecting line
    mark_inset(ax, axins, loc1=2, loc2=4, fc="none", ec="0.5")

for T, U, V in iterate_containers(Ts, Us, Vs):
    name = f"T={T}/U={U}/V={V}"
    data, data_real, w_lin, res = cf.resolvent_data(f"{folder}{name}", name_suffix, plot_lower_lim, plot_upper_lim, number_of_values=20000, xp_basis=use_XP)
    ax.plot(w_lin, data, label=name_suffix)
    if createZoom:
        axins.plot(w_lin, data)

name_suffix = "AFM"
for T, U, V in iterate_containers(Ts, Us, Vs):
    name = f"T={T}/U={U}/V={V}"
    data, data_real, w_lin, res = cf.resolvent_data(f"{folder}{name}", name_suffix, plot_lower_lim, plot_upper_lim, number_of_values=20000, xp_basis=use_XP)
    ax.plot(w_lin, data, linestyle="--", label=name_suffix)
    if createZoom:
        axins.plot(w_lin, data)

res.mark_continuum(ax)
if createZoom:
    res.mark_continuum(axins)
legend = ax.legend()

ax.set_xlim(plot_lower_lim, plot_upper_lim)
ax.set_xlabel(r"$z / t$")
ax.set_ylabel(r"Spectral density / a.u.")
fig.tight_layout()

import os
plt.savefig(f"python/build/{os.path.basename(__file__).split('.')[0]}.pdf")
plt.show()
