import numpy as np
import matplotlib.pyplot as plt
import lib.continued_fraction as cf
import lib.plot_settings as ps

use_XP = True
fig, ax = plt.subplots()

plotter = ps.CURVEFAMILY(3, axis=ax)
plotter.set_individual_colors("nice")
plotter.set_individual_linestyles(["-", "-.", "--"])

plot_lower_lim = -12
plot_upper_lim = 12

name_suffix = "phase_SC"

for val, disc in zip([900, 3000, 6000], ["900", "3k", "6k"]):
    folder = f"data/modes/cube/dos_{disc}/T=0.0/U=0.0/V=0.0"
    data, data_real, w_lin, res = cf.resolvent_data(folder, name_suffix, plot_lower_lim, plot_upper_lim, 
                                                    number_of_values=2000, xp_basis=use_XP, imaginary_offset=0)
    plotter.plot(w_lin, data, label=f"$N_\\gamma = {val}$")

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
