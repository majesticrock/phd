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
U = -2.5
V = -0.1

folders = ["../data/modes/square/dos_3k/", "../data/modes/cube/dos_3k/"]
labels = ["Square", "Simple cubic"]
name_suffix = "phase_SC"
fig, axs = plt.subplots(nrows=2, figsize=(6.4, 4.8), sharex=True, gridspec_kw=dict(hspace=0, wspace=0))

plotters = np.empty(2, dtype=ps.CURVEFAMILY)

for i in range(2):
    plotters[i] = ps.CURVEFAMILY(4, axis=axs[i])
    plotters[i].set_individual_colors("nice2")
    plotters[i].set_individual_dashes()

    name = f"T={T}/U={U}/V={V}"
    data_imag, data, w_log, res = cf.resolvent_data_log_z(f"{folders[i]}{name}", name_suffix, lower_edge=0, begin_offset=3e-3, range=0.2, xp_basis=True, imaginary_offset=0, number_of_values=1000, messages=False)
    
    data = np.log(data)
    plotters[i].plot(w_log, data, label=f"Data")
    popt, pcov = curve_fit(func_ln, w_log, data)
    plotters[i].plot(w_log, func_ln(w_log, *popt), label=f"Fit")
    axs[i].text(0.05, 0.4, f"$a={popt[0]:.4f}\pm{np.sqrt(pcov[0][0]):.4f}$", transform = axs[i].transAxes)
    axs[i].text(0.05, 0.3, f"$b={popt[1]:.4f}\pm{np.sqrt(pcov[1][1]):.4f}$", transform = axs[i].transAxes)

axs[0].legend()
axs[1].set_xlabel(r"$\ln(\omega / t)$")
axs[1].text(-0.08, 1, r"$\ln(\Re[\mathcal{G}_\mathrm{Phase}](\omega))$", va='center', ha='center', rotation='vertical', transform = axs[1].transAxes)
fig.tight_layout()

import os
plt.savefig(f"plots/{os.path.basename(__file__).split('.')[0]}.pdf")
#plt.show()