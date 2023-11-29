import numpy as np
import matplotlib.pyplot as plt

import os, sys
if os.name == "nt":
    sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + r"\python")
else:
    sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + "/python")
import lib.continued_fraction as cf
import lib.plot_settings as ps

from scipy.optimize import curve_fit
def func_ln(x, a, b):
    return a * x + b

T = 0.
Us = [-2., -2.5]
V = -0.1
use_xp = True

folders = ["data/modes/square/dos_3k/", "data/modes/cube/dos_3k/"]
labels = ["Square", "Simple cubic"]
name_suffix = "phase_SC"
fig, axs = plt.subplots(nrows=2, figsize=(6.4, 4.8), sharex=True)

plotters = np.empty(2, dtype=ps.CURVEFAMILY)
plot_lower_lim = 0.005
plot_upper_lim = plot_lower_lim + 0.2

for i in range(2):
    plotters[i] = ps.CURVEFAMILY(4, axis=axs[0])
    plotters[i].set_individual_colors("nice2")
    plotters[i].set_individual_linestyles(["-", "--"])

    name = f"T={T}/U={Us[i]}/V={V}"
    data_imag, data, w_lin, res = cf.resolvent_data(f"{folders[i]}{name}", name_suffix, plot_lower_lim, plot_upper_lim, xp_basis=use_xp, imaginary_offset=0)

    w_log = np.log(w_lin.real)
    plotters[i].plot(w_log, np.log(data), label=f"{labels[i]} data")
    popt, pcov = curve_fit(func_ln, w_log, np.log(data))
    plotters[i].plot(w_log, func_ln(w_log, *popt), label=f"{labels[i]} fit")
    axs[i].legend()
    axs[i].text(0.05, 0.35, f"$a={popt[0]:.4f}\pm{np.sqrt(pcov[0][0]):.4f}$")
    axs[i].text(0.05, 0.3, f"$b={popt[1]:.4f}\pm{np.sqrt(pcov[1][1]):.4f}$")
    
axs[0].set_xlabel(r"$\ln(\omega)$")
axs[0].set_ylabel(r"$\ln(\Re G^\mathrm{ret}(\omega))$")
axs[1].set_ylabel(r"$\ln(\Re G^\mathrm{ret}(\omega))$")


fig.tight_layout()

import os
#plt.savefig(f"plots/{os.path.basename(__file__).split('.')[0]}.pdf")
plt.show()