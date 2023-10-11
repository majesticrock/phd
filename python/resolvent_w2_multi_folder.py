import numpy as np
import matplotlib.pyplot as plt
import lib.continued_fraction as cf
from lib.iterate_containers import iterate_containers
# Calculates the resolvent in w^2

prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']

Ts = np.array([0.])
Us = np.array([-2.0])
Vs = np.array([0.25, 0.1, 0.05, 0.01, 0.005, 0.001, 0])

folder = "data/modes/square/dos_900/"
element_names = ["a", "a+b", "a+ib"]
fig, ax = plt.subplots()

#cax.set_xscale("log")
#ax.set_yscale("log")

plot_upper_lim = 8.5
name_suffix = "phase_SC"

peak_positions = np.zeros(len(Vs))
counter = 0
for T, U, V in iterate_containers(Ts, Us, Vs):
    name = f"T={T}/U={U}_V={V}"
    data, data_real, w_lin, res = cf.resolvent_data(f"{folder}{name}", name_suffix, 1e-4, plot_upper_lim, number_of_values=30000, imaginary_offset=1e-6)
    #ax.plot(w_lin, data_real, linewidth=(plt.rcParams["lines.linewidth"]), linestyle="--", label=f"$V={V}$")
    
    peak_positions[counter] = w_lin[np.argmax(data)]
    counter += 1

ax.plot(Vs, peak_positions)
#res.mark_continuum(ax)
#legend = plt.legend(loc=8)

ax.set_xlabel(r"$z / t$")
ax.set_ylabel(r"Spectral density / a.u.")
fig.tight_layout()

import os
plt.savefig(f"python/build/{os.path.basename(__file__).split('.')[0]}.pdf")
plt.show()
