import numpy as np
import matplotlib.pyplot as plt
import lib.continued_fraction as cf
import lib.plot_settings as ps
from lib.iterate_containers import *
# Calculates the resolvent in w^2

prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']

realPart = False
both = False

Ts = np.array([0.])
Us = np.array([-2.0])
Vs = np.array(["0.002", "0.001", "0.0005", "0.0001", "0.00005", "0.0"])#"0.25", "0.2", "0.15", "0.1", "0.05", "0.04", "0.03", "0.02", "0.01", "0.007", "0.005", "0.004", 

folder = "data/modes/square/dos_900/"
element_names = ["a", "a+b", "a+ib"]
fig, ax = plt.subplots()

if realPart or both:
    ax.set_xscale("log")
    ax.set_yscale("symlog")
else:
    ax.set_yscale("log")

plot_upper_lim = 9
name_suffix = "phase_SC"

realPlotter = ps.CURVEFAMILY(total_size(Ts, Us, Vs), axis=ax, allow_cycle=True)
realPlotter.set_individual_colors("nice")

plotter = ps.CURVEFAMILY(total_size(Ts, Us, Vs), axis=ax, allow_cycle=True)
plotter.set_individual_colors("nice")
plotter.set_shared_linestyle("--")

for T, U, V in iterate_containers(Ts, Us, Vs):
    name = f"T={T}/U={U}_V={V}"
    data, data_real, w_lin, res = cf.resolvent_data(f"{folder}{name}", name_suffix, 1e-6, plot_upper_lim, number_of_values=20000, imaginary_offset=1e-6)
    if realPart or both:
        realPlotter.plot(w_lin, data_real, linewidth=(plt.rcParams["lines.linewidth"]), label=f"$V={V}$")
    if not realPart or both:
        plotter.plot(w_lin, data, linewidth=(plt.rcParams["lines.linewidth"]), label=f"$V={V}$")

#res.mark_continuum(ax)
#legend = plt.legend(loc=8)

ax.set_xlabel(r"$z / t$")
ax.set_ylabel(r"Spectral density / a.u.")
#fig.tight_layout()

import os
plt.savefig(f"python/build/{os.path.basename(__file__).split('.')[0]}.pdf")
plt.show()
