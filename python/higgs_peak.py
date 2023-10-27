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

folder = "data/modes/square/dos_900/"
fig, ax = plt.subplots()

#ax.set_xscale("log")
#ax.set_yscale("log")

name = f"T={T}/U={U}/V={V}"
phase_data, data_real, w_lin, res = cf.resolvent_data_log_z(f"{folder}{name}", "phase_SC", begin_offset=1e-5, range=0.1, number_of_values=20000, xp_basis=use_XP)
higgs_data, data_real, w_lin, res = cf.resolvent_data_log_z(f"{folder}{name}", "higgs_SC", begin_offset=1e-5, range=0.1, number_of_values=20000, xp_basis=use_XP)

diff_data = higgs_data - phase_data
ax.plot(w_lin, higgs_data, linewidth=(plt.rcParams["lines.linewidth"]), linestyle="-", label="Higgs")
#ax.plot(w_lin, diff_data, linewidth=(plt.rcParams["lines.linewidth"]), linestyle="--", label="Higgs - Phase")
#ax.plot(w_lin, -data_real, label="Real part")

from scipy.optimize import curve_fit
def func(x, a, b):
    return a * x + b

def func2(x, a, b, c):
    return a * (np.tanh(b*x - c) + 1)

cut = 5000
popt, pcov = curve_fit(func, w_lin[:cut], diff_data[:cut])
ax.plot(w_lin, func(w_lin, *popt), "k--", label="Fit")

#popt, pcov = curve_fit(func2, w_lin, -data_real, p0=(11, -1, 10))
#ax.plot(w_lin, func2(w_lin, *popt), "k:", label="Fit")

legend = plt.legend()
ax.add_artist(legend)

ax.set_xlabel(r"$\ln(z - z_-) / t$")
ax.set_ylabel(r"Spectral density / a.u.")
fig.tight_layout()

import os
plt.savefig(f"python/build/{os.path.basename(__file__).split('.')[0]}.pdf")
plt.show()
