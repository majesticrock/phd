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
    plotters[i].set_individual_colors(["blue", "black", "orange"])
    plotters[i].set_individual_dashes([(1, 0), (1, 1), (5, 5)])

    name = f"T={T}/U={U}/V={V}"
    phase_imag, phase_real, w_log, res = cf.resolvent_data_log_z(f"{folders[i]}{name}", "phase_SC", range=0.05, begin_offset=1e-3, number_of_values=10000, imaginary_offset=0, xp_basis=True, messages=False)
    higgs_imag, higgs_real, w_log, res = cf.resolvent_data_log_z(f"{folders[i]}{name}", "higgs_SC", range=0.05, begin_offset=1e-3, number_of_values=10000, imaginary_offset=0, xp_basis=True, messages=False)

    diff_data = np.log(higgs_imag - phase_imag)
    plotters[i].plot(w_log, diff_data, label="Higgs - Phase", linewidth=1.25*plt.rcParams["lines.linewidth"])
    plotters[i].plot(w_log, np.log(higgs_imag), label="Higgs", linewidth=1.75*plt.rcParams["lines.linewidth"])

    from scipy.optimize import curve_fit
    def func(x, a, b):
        return a * x + b

    popt, pcov = curve_fit(func, w_log, diff_data)
    axs[i].text(0.05, 0.4, f"$a={popt[0]:.4f}$", transform = axs[i].transAxes)
    axs[i].text(0.05, 0.3, f"$b={popt[1]:.4f}$", transform = axs[i].transAxes)
    plotters[i].plot(w_log, func(w_log, *popt), label="Fit")

axs[0].legend(loc="upper right")

axs[1].set_xlabel(r"$\ln((\omega - \omega_0) / t)$")
axs[1].text(-0.08, 1, r"$\ln(-\Im \mathcal{G}_\mathrm{Higgs}(\omega))$", va='center', ha='center', rotation='vertical', transform = axs[1].transAxes)
fig.tight_layout()

plt.savefig(f"plots/{os.path.basename(__file__).split('.')[0]}.pdf")
#plt.show()
