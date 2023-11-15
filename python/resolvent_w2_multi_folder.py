import numpy as np
import matplotlib.pyplot as plt
import lib.continued_fraction as cf
import lib.plot_settings as ps
from lib.iterate_containers import *
from lib.extract_key import * 
# Calculates the resolvent in w^2

prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']

realPart = False
both = False

Ts = np.array([0.0])
Us = np.array([-2.0, -2.5, -3., -3.5])
Vs = np.array([-0.1])

folder = "data/modes/cube/dos_6k/"
fig, ax = plt.subplots()

if realPart or both:
    ax.set_xscale("log")
    ax.set_yscale("symlog")
else:
    ax.set_yscale("log")

plot_upper_lim = 1
name_suffix = "phase_SC"

realPlotter = ps.CURVEFAMILY(total_size(Ts, Us, Vs), axis=ax, allow_cycle=True)
realPlotter.set_individual_colors("default")

plotter = ps.CURVEFAMILY(total_size(Ts, Us, Vs), axis=ax, allow_cycle=True)
plotter.set_individual_colors("default")
plotter.set_individual_linestyles()
#plotter.set_shared_linestyle("-")

for T, U, V in iterate_containers(Ts, Us, Vs):
    name = f"T={T}/U={U}/V={V}"
    if len(Ts) > 1:
        label = f"$T={T}$"
    elif len(Us) > 1:
        label = f"$U={U}$"
    elif len(Vs) > 1:
        label = f"$V={V}$"
    
    data, data_real, w_lin, res = cf.resolvent_data(f"{folder}{name}", name_suffix, 0, upper_edge=plot_upper_lim, number_of_values=20000, imaginary_offset=1e-6, xp_basis=True)
    print(w_lin[np.argmax(data)])
    if realPart or both:
        realPlotter.plot(w_lin, data_real, label=label)
    if not realPart or both:
        plotter.plot(w_lin, data, label=label)

#res.mark_continuum(ax)
legend = plt.legend(loc=8)

ax.set_xlabel(r"$z / t$")
ax.set_ylabel(r"Spectral density / a.u.")
#fig.tight_layout()

import os
plt.savefig(f"python/build/{os.path.basename(__file__).split('.')[0]}.pdf")
plt.show()
