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
V = 0.0

folders = ["../data/modes/square/dos_3k/", "../data/modes/cube/dos_3k/"]

fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(12.8, 4.8), sharex=True, sharey=True, gridspec_kw=dict(hspace=0, wspace=0))

for j, V, g_name, name_suffix in [(0, -0.1, "Phase", "phase_SC"), (1, 0.0, "Higgs", "higgs_SC")]:
    plotters = np.empty(2, dtype=ps.CURVEFAMILY)
    for i in range(2):
        plotters[i] = ps.CURVEFAMILY(4, axis=axs[j][i])
        plotters[i].set_individual_colors("nice2")
        plotters[i].set_individual_dashes()

        name = f"T={T}/U={U}/V={V}"
        data_imag, data, w_log, res = cf.resolvent_data_log_z(f"{folders[i]}{name}", name_suffix, lower_edge=0, begin_offset=3e-3, range=0.2, xp_basis=True, imaginary_offset=0, number_of_values=1000, messages=False)

        data = np.log(data)
        plotters[i].plot(w_log, data, label=f"Data")
        popt, pcov = curve_fit(func_ln, w_log, data)
        plotters[i].plot(w_log, func_ln(w_log, *popt), label=f"Fit")
        axs[j][i].text(0.05, 0.4, f"$a={popt[0]:.4f}\pm{np.sqrt(pcov[0][0]):.4f}$", transform = axs[j][i].transAxes)
        axs[j][i].text(0.05, 0.29, f"$b={popt[1]:.4f}\pm{np.sqrt(pcov[1][1]):.4f}$", transform = axs[j][i].transAxes)
    axs[j][0].set_ylabel(rf"$\ln(\Re[\mathcal{{G}}_\mathrm{{{g_name}}}](\omega))$")

axs[0][0].title.set_text("Square")
axs[0][1].title.set_text("Simple cubic")
axs[0][0].legend()
axs[1][0].set_xlabel(r"$\ln(\omega / t)$")
axs[1][1].set_xlabel(r"$\ln(\omega / t)$")

axs[0][1].text(-2.8, 8.5, "(a) SC\n$V=-0.1t$")
axs[1][1].text(-2.8, 8.5, "(b) SC/CDW\n$V=0$")

fig.tight_layout()

plt.savefig(f"plots/{os.path.basename(__file__).split('.')[0]}.pdf")
#plt.show()