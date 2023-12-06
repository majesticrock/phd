import numpy as np
import matplotlib.pyplot as plt

import os, sys
if os.name == "nt":
    sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + r"\python")
else:
    sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + "/python")
import lib.continued_fraction as cf
import lib.plot_settings as ps
import lib.resolvent_peak as rp

T = 0.
U = -2.5
V = -0.1

folders = ["../data/modes/square/dos_3k/", "../data/modes/cube/dos_3k/"]
labels = ["Square", "Simple cubic"]
name_suffix = "CDW"
fig, axs = plt.subplots(nrows=2, figsize=(6.4, 4.8), sharex=True, gridspec_kw=dict(hspace=0, wspace=0))

plotters = np.empty(2, dtype=ps.CURVEFAMILY)
plot_lower_lim = 0.006
plot_upper_lim = plot_lower_lim + 0.1

for i in range(2):
    plotters[i] = ps.CURVEFAMILY(4, axis=axs[i])
    plotters[i].set_individual_colors("nice2")
    plotters[i].set_individual_dashes()

    name = f"T={T}/U={U}/V={V}"
    peak = rp.Peak(f"{folders[i]}{name}", name_suffix, initial_search_bounds=(0.3, 0.8))
    #print(peak.peak_position)
    peak.improved_peak_position(x0_offset=1e-1, gradient_epsilon=1e-12)
    popt, pcov, w_space, y_data = peak.fit_real_part(range=0.05, begin_offset=1e-10, reversed=False)
    
    axs[i].text(0.05, 0.4, f"$a={popt[0]:.4f}\pm{np.sqrt(pcov[0][0]):.4f}$", transform = axs[i].transAxes)
    axs[i].text(0.05, 0.3, f"$b={popt[1]:.4f}\pm{np.sqrt(pcov[1][1]):.4f}$", transform = axs[i].transAxes)
    plotters[i].plot(w_space, y_data, label=f"{labels[i]} data")
    plotters[i].plot(w_space, rp.linear_function(w_space, *popt), "k--", label=f"{labels[i]} fit")
    axs[i].set_ylim(0, 22.5)#
    axs[i].legend()
    
axs[1].set_xlabel(r"$\ln((\omega - \omega_0) / t)$")
axs[1].text(-0.08, 1, r"$\ln(\Re[\mathcal{G}_\mathrm{CDW}](\omega - \omega_0))$", va='center', ha='center', rotation='vertical', transform = axs[1].transAxes)
fig.tight_layout()

import os
plt.savefig(f"plots/{os.path.basename(__file__).split('.')[0]}.pdf")
#plt.show()