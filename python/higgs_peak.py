import numpy as np
import matplotlib.pyplot as plt
import lib.continued_fraction as cf
from lib.iterate_containers import iterate_containers
# Calculates the resolvent in w^2

prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']

T = 0.
U = -2.0
V = -0.1

use_XP = True

folder = "data/modes/square/dos_2500/"
fig, ax = plt.subplots()

#ax.set_xscale("log")
#ax.set_yscale("log")

name = f"T={T}/U={U}/V={V}"
phase_data, data_real, w_lin, res = cf.resolvent_data_log_z(f"{folder}{name}", "phase_SC", range=0.45, begin_offset=1e-4, number_of_values=20000, xp_basis=use_XP)
higgs_data, data_real, w_lin, res = cf.resolvent_data_log_z(f"{folder}{name}", "higgs_SC", range=0.45, begin_offset=1e-4, number_of_values=20000, xp_basis=use_XP)

diff_data = np.log(higgs_data )
ax.plot(w_lin, diff_data, linewidth=(plt.rcParams["lines.linewidth"]), linestyle="-", label="Higgs")
#ax.plot(w_lin, diff_data, linewidth=(plt.rcParams["lines.linewidth"]), linestyle="--", label="Higgs - Phase")
#ax.plot(w_lin, -data_real, label="Real part")

from scipy.optimize import curve_fit
def func(x, a, b):
    return a * x + b

popt, pcov = curve_fit(func, w_lin, diff_data)
print(popt[0], " +/- ", np.sqrt(pcov[0][0]))
print(popt[1], " +/- ", np.sqrt(pcov[1][1]))
ax.plot(w_lin, func(w_lin, *popt), "k--", label="Fit")

legend = plt.legend()
ax.add_artist(legend)

ax.set_xlabel(r"$\ln(z - z_-) / t$")
ax.set_ylabel(r"Spectral density / a.u.")
fig.tight_layout()

import os
plt.savefig(f"python/build/{os.path.basename(__file__).split('.')[0]}.pdf")
plt.show()
