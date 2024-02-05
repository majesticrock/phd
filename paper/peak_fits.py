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

folders = ["../data/modes/square/dos_3k/", "../data/modes/cube/dos_3k/"]

fig, axs = plt.subplots(nrows=4, ncols=2, figsize=(12.8, 9.6), sharey="row", sharex=True, gridspec_kw=dict(hspace=0, wspace=0))
def plot_general_peak(plot_axs, T, U, V, name_suffix, initial_search_bounds, text="", reversed=False):
    plotters = np.empty(2, dtype=ps.CURVEFAMILY)
    
    for i in range(2):
        plotters[i] = ps.CURVEFAMILY(2, axis=plot_axs[i])
        plotters[i].set_individual_colors("nice2")
        plotters[i].set_individual_dashes()

        name = f"T={T}/U={U[i] if hasattr(U, '__len__') else U}/V={V}"
        peak = rp.Peak(f"{folders[i]}{name}", name_suffix, initial_search_bounds=initial_search_bounds)

        peak.improved_peak_position(xtol=1e-13)
        popt, pcov, w_space, y_data = peak.fit_real_part(range=0.05, begin_offset=1e-10, reversed=reversed)

        plot_axs[i].text(0.05, 0.4, f"$a={popt[0]:.4f}$", transform = plot_axs[i].transAxes)
        plot_axs[i].text(0.05, 0.3, f"$b={popt[1]:.4f}$", transform = plot_axs[i].transAxes)
        plotters[i].plot(w_space, y_data, label="Data")
        plotters[i].plot(w_space, rp.linear_function(w_space, *popt), label="Fit")

        plot_axs[i].text(0.5, 0.7, f"{text}$U={U[i] if hasattr(U, '__len__') else U}, V={V}$", transform = plot_axs[i].transAxes)

# CDW in SC
plot_general_peak(axs[0], T=0.0, U=-2.5, V=-0.1, name_suffix="CDW", initial_search_bounds=(0.3, 0.8), reversed=False, text="SC: ")
axs[0][0].set_ylabel(r"$\ln(\Re[\mathcal{G}_\mathrm{CDW}](\omega) \cdot t)$")
# SC in CDW
plot_general_peak(axs[1], T=0.0, U=-2.5, V=0.1, name_suffix="higgs_SC", initial_search_bounds=(0.5, 1.5), reversed=False, text="CDW: ")
axs[1][0].set_ylabel(r"$\ln(\Re[\mathcal{G}_\mathrm{Higgs}](\omega) \cdot t)$")
#AFM in CDW
plot_general_peak(axs[2], T=0.0, U=[3.9, 5.9], V=1.0, name_suffix="AFM", initial_search_bounds=(2., "lower_edge"), reversed=True, text="CDW: ")
axs[2][0].set_ylabel(r"$\ln(\Re[\mathcal{G}_\mathrm{AFM}](\omega) \cdot t)$")
#AFM in AFM
plot_general_peak(axs[3], T=0.0, U=[4.1, 6.1], V=1.0, name_suffix="AFM", initial_search_bounds=(2., "lower_edge"), reversed=True, text="AFM: ")
axs[3][0].set_ylabel(r"$\ln(\Re[\mathcal{G}_\mathrm{AFM}](\omega) \cdot t)$")


legend = axs[0][1].legend(loc='upper center', bbox_to_anchor=(0., 1.25), ncol=2, shadow=True)
axs[0][0].title.set_text("Square")
axs[0][1].title.set_text("Simple cubic")
axs[len(axs)-1][0].set_xlabel(r"$\ln((\omega - \omega_0) / t)$")
axs[len(axs)-1][1].set_xlabel(r"$\ln((\omega - \omega_0) / t)$")
fig.tight_layout()

plt.savefig(f"plots/{os.path.basename(__file__).split('.')[0]}.pdf")
#plt.show()