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

folders = ["../data/modes/square/dos_6000", "../data/modes/cube/dos_6000"]
folder_extensions = ["_SC/", "_CDW/", "/"]
name_suffices = ["higgs_SC", "higgs_CDW", "higgs_SC"]

nrows = 3
fig, axs = plt.subplots(nrows=3, ncols=2, figsize=(12.8, nrows * 0.5 * 4.8), sharex=True, sharey=True, gridspec_kw=dict(hspace=0, wspace=0))

plotters = np.empty((nrows, 2), dtype=ps.CURVEFAMILY)
for i in range(2):
    for j in range(nrows):
        plotters[j][i] = ps.CURVEFAMILY(3, axis=axs[j][i])
        plotters[j][i].set_individual_colors(["blue", "black", "orange"])
        plotters[j][i].set_individual_dashes([(1, 0), (1, 1), (5, 5)])

        name = f"T={T}/U={U}/V={V}"
        phase_imag, phase_real, w_log, res = cf.resolvent_data_log_z(f"{folders[i]}{folder_extensions[j]}{name}", "phase_SC", range=0.08, begin_offset=1e-4, number_of_values=10000, imaginary_offset=0, xp_basis=True, messages=False)
        higgs_imag, higgs_real, w_log, res = cf.resolvent_data_log_z(f"{folders[i]}{folder_extensions[j]}{name}", name_suffices[j], range=0.08, begin_offset=1e-4, number_of_values=10000, imaginary_offset=0, xp_basis=True, messages=False)

        diff_data = np.log(higgs_imag - phase_imag)
        plotters[j][i].plot(w_log, diff_data, label=r"$\mathcal{A}_i (\omega) - \mathcal{A}_\mathrm{Phase} (\omega)$", linewidth=1.25*plt.rcParams["lines.linewidth"])
        plotters[j][i].plot(w_log, np.log(higgs_imag), label=r"$\mathcal{A}_i (\omega)$", linewidth=1.75*plt.rcParams["lines.linewidth"])

        from scipy.optimize import curve_fit
        def func(x, a, b):
            return a * x + b

        popt, pcov = curve_fit(func, w_log, diff_data)
        axs[j][i].text(0.05, 0.2, f"$a={popt[0]:.4f}$", transform = axs[j][i].transAxes)
        axs[j][i].text(0.05, 0.1, f"$b={popt[1]:.4f}$", transform = axs[j][i].transAxes)
        plotters[j][i].plot(w_log, func(w_log, *popt), label="Fit")

axs[0][0].text(-4.6, 1.15, "(a.1)\n$\Delta_\\mathrm{CDW} = 0$")
axs[1][0].text(-4.6, 1.15, "(b.1)\n$\Delta_\\mathrm{SC} = 0$")
axs[2][0].text(-4.6, 1.15, "(c.1)\n$\Delta_\\mathrm{SC} = \Delta_\\mathrm{CDW}$")
axs[0][1].text(-4.6, 1.15, "(a.2)\n$\Delta_\\mathrm{CDW} = 0$")
axs[1][1].text(-4.6, 1.15, "(b.2)\n$\Delta_\\mathrm{SC} = 0$")
axs[2][1].text(-4.6, 1.15, "(c.2)\n$\Delta_\\mathrm{SC} = \Delta_\\mathrm{CDW}$")

axs[0][0].title.set_text("Square")
axs[0][1].title.set_text("Simple cubic")
legend = axs[0][1].legend(loc='upper center', bbox_to_anchor=(0., 1.3))

axs[nrows-1][0].set_xlabel(r"$\ln((\omega - 2 \Delta_\mathrm{tot}) / t)$")
axs[nrows-1][1].set_xlabel(r"$\ln((\omega - 2 \Delta_\mathrm{tot}) / t)$")

axs[0][0].set_ylabel(r"$\ln(\mathcal{A}_\mathrm{Higgs}(\omega) \cdot t)$")
axs[1][0].set_ylabel(r"$\ln(\mathcal{A}_\mathrm{CDW}  (\omega) \cdot t)$")
axs[2][0].set_ylabel(r"$\ln(\mathcal{A}_\mathrm{Higgs}(\omega) \cdot t)$")
fig.tight_layout()

plt.savefig(f"plots/{os.path.basename(__file__).split('.')[0]}.pdf")
#plt.show()
