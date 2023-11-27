import numpy as np
import matplotlib.pyplot as plt
import lib.continued_fraction as cf
from lib.iterate_containers import naming_scheme
import lib.plot_settings as ps
# Calculates the resolvent in w^2

prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']

Ts = np.array([0.])
Us = np.array([-2.0])
Vs = np.array([-0.1])

use_XP = True

folder = "data/modes/cube/dos_3k/"
fig, ax = plt.subplots()

#ax.set_xscale("log")
#ax.set_yscale("log")
#ax.set_ylim(0, 1)

plotter = ps.CURVEFAMILY(3, axis=ax)
plotter.set_individual_colors("nice")
plotter.set_individual_linestyles(["-", "-.", "--"])

plot_lower_lim = 0
plot_upper_lim = 12

name_suffix = "phase_SC"
for name in naming_scheme(Ts, Us, Vs):
    data, data_real, w_lin, res = cf.resolvent_data(f"{folder}{name}", name_suffix, plot_lower_lim, plot_upper_lim, 
                                                    number_of_values=20000, xp_basis=use_XP, imaginary_offset=1e-6)
    plotter.plot(w_lin, data, label="Phase")

name_suffix = "higgs_SC"
for name in naming_scheme(Ts, Us, Vs):
    data, data_real, w_lin, res = cf.resolvent_data(f"{folder}{name}", name_suffix, plot_lower_lim, plot_upper_lim, 
                                                    number_of_values=20000, xp_basis=use_XP, imaginary_offset=1e-6)
    plotter.plot(w_lin, data, label="Higgs")
    
name_suffix = "CDW"
for name in naming_scheme(Ts, Us, Vs):
    data, data_real, w_lin, res = cf.resolvent_data(f"{folder}{name}", name_suffix, plot_lower_lim, plot_upper_lim, 
                                                    number_of_values=20000, xp_basis=use_XP, imaginary_offset=1e-6)
    plotter.plot(w_lin, data, label=name_suffix)

res.mark_continuum(ax)
legend = plt.legend()
ax.add_artist(legend)

ax.set_xlim(plot_lower_lim, plot_upper_lim)
ax.set_xlabel(r"$z / t$")
ax.set_ylabel(r"Spectral density / a.u.")
fig.tight_layout()

import os
plt.savefig(f"python/build/{os.path.basename(__file__).split('.')[0]}.pdf")
plt.show()
