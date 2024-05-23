import numpy as np
import matplotlib.pyplot as plt

import os, sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import lib.continued_fraction as cf
import lib.plot_settings as ps
# Calculates the resolvent in w^2

prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']

use_XP = True

folder = "data/continuum/test"
fig, ax = plt.subplots()

ax.set_yscale("log")
#ax.set_yscale("symlog")
#ax.set_ylim(-0.05, 100.)

plotter = ps.CURVEFAMILY(6, axis=ax)
plotter.set_individual_colors("nice")
plotter.set_individual_linestyles(["-", "-.", "--", "-", "--", ":"])
#plotter.set_individual_dashes()

plot_lower_lim = 0
plot_upper_lim = 0.07

name_suffix = "phase_SC"
data, data_real, w_lin, res = cf.resolvent_data(f"{folder}", name_suffix, plot_lower_lim, plot_upper_lim, 
                                                    number_of_values=20000, xp_basis=use_XP, imaginary_offset=1e-6, ingore_first=5)
plotter.plot(w_lin, data, label="Phase")

name_suffix = "higgs_SC"
data, data_real, w_lin, res = cf.resolvent_data(f"{folder}", name_suffix, plot_lower_lim, plot_upper_lim, 
                                                    number_of_values=20000, xp_basis=use_XP, imaginary_offset=1e-6, ingore_first=5)
plotter.plot(w_lin, data, label="Higgs")

res.mark_continuum(ax)
legend = plt.legend()
ax.add_artist(legend)

ax.set_xlim(plot_lower_lim, plot_upper_lim)
ax.set_xlabel(r"$\omega [\mathrm{meV}]$")
ax.set_ylabel(r"$\mathcal{A} (\omega) [1 / \mathrm{meV}]$")
fig.tight_layout()

import os
plt.savefig(f"python/build/{os.path.basename(__file__).split('.')[0]}.pdf")
plt.show()
