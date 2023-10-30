import numpy as np
import matplotlib.pyplot as plt
import lib.continued_fraction as cf
from lib.iterate_containers import iterate_containers
# Calculates the resolvent in w^2

prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']

Ts = np.array([0.0])
Us = np.array([-10.0])
Vs = np.array([-0.1])

use_XP = True

folder = "data/modes/square/dos_2500/"
name_suffix = "phase_SC"
element_names = ["a", "a+b", "a+ib"]
fig, ax = plt.subplots()

#ax.set_xscale("log")
ax.set_yscale("log")

plot_lower_lim = 0
plot_upper_lim = 8

for T, U, V in iterate_containers(Ts, Us, Vs):
    name = f"T={T}/U={U}/V={V}"
    data, data_real, w_lin, res = cf.resolvent_data(f"{folder}{name}", name_suffix, plot_lower_lim, plot_upper_lim, number_of_values=20000, xp_basis=use_XP)
    ax.plot(w_lin, data, linewidth=(plt.rcParams["lines.linewidth"]), label=name_suffix)

name_suffix = "higgs_SC"
for T, U, V in iterate_containers(Ts, Us, Vs):
    name = f"T={T}/U={U}/V={V}"
    data, data_real, w_lin, res = cf.resolvent_data(f"{folder}{name}", name_suffix, plot_lower_lim, plot_upper_lim, number_of_values=20000, xp_basis=use_XP)
    ax.plot(w_lin, data, linewidth=(plt.rcParams["lines.linewidth"]), linestyle="--", label=name_suffix)
    #ax.plot(w_lin, 0.75*0.1*(np.log(w_lin - np.sqrt(res.roots[0])))**2 )

res.mark_continuum(ax)
legend = plt.legend()

#import matplotlib.lines as mlines
#dummy_lines = []
#dummy_lines.append(mlines.Line2D([],[], color="k", linestyle="-"))
#dummy_lines.append(mlines.Line2D([],[], color="k", linestyle="--", linewidth=2*plt.rcParams["lines.linewidth"]))
#legend_extra = plt.legend([dummy_lines[i] for i in [0,1]], [r"Amplitude", r"Phase"], loc="upper center")

ax.add_artist(legend)
#ax.add_artist(legend_extra)

ax.set_xlabel(r"$z / t$")
ax.set_ylabel(r"Spectral density / a.u.")
fig.tight_layout()

import os
plt.savefig(f"python/build/{os.path.basename(__file__).split('.')[0]}.pdf")
plt.show()
