import numpy as np
import matplotlib.pyplot as plt

import os, sys
if os.name == "nt":
    sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + r"\python")
else:
    sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + "/python")

import lib.continued_fraction as cf
from lib.iterate_containers import naming_scheme_tuples
import lib.plot_settings as ps
from lib.create_zoom import *

T = 0.0
U = -2.5
V = -0.1

folders = ["../data/modes/square/", "../data/modes/cube/"]
folder_extensions = ["dos_300/", "dos_1500/", "dos_3k/", "dos_6k/"]
nrows = 4
ncols = 2
# ax = axs[row][col]
fig, axs = plt.subplots(nrows=nrows, ncols=ncols, figsize=(12.8, 8), sharey=True, sharex="col", gridspec_kw=dict(hspace=0, wspace=0))

plotters = np.empty((nrows, ncols), dtype=ps.CURVEFAMILY)
for i in range(nrows):
    for j in range(ncols):
        axs[i][j].set_ylim(0, .99)
        plotters[i][j] = ps.CURVEFAMILY(4, axis=axs[i][j])
        plotters[i][j].set_individual_colors("nice")
        plotters[i][j].set_individual_dashes([ [1,0], [3,6,3,0], [1.5, 2], [3,5,1.5,5,3,0] ])

plot_lower_lim = -0.05
plot_upper_lim = 4.25

name_suffices = ["phase_SC", "higgs_SC", "CDW", "AFM"]
labels = ["Phase", "Higgs", "CDW", "AFM"]

for j, folder in enumerate(folders):
    usage_upper_lim = 2 * plot_upper_lim if j == 0 else 3 * plot_upper_lim
    name = f"T={T}/U={U}/V={V}"
    for i, folder_extension in enumerate(folder_extensions):
        resolvents = np.empty(len(name_suffices), dtype=cf.ContinuedFraction)
        for k, (name_suffix, label) in enumerate(zip(name_suffices, labels)):
            data, data_real, w_lin, resolvents[k] = cf.resolvent_data(f"{folder}{folder_extension}{name}", name_suffix, plot_lower_lim, usage_upper_lim, 
                                                            number_of_values=10000, xp_basis=True, imaginary_offset=1e-5, messages=False, ingore_first=10*i+5)
            plotters[i][j].plot_with_peak(w_lin, data, label=r"$\mathcal{A}_\mathrm{" + label + r"} (\omega)$")
            
        axs[i][j].set_xlim(plot_lower_lim, usage_upper_lim)
        resolvents[0].mark_continuum(axs[i][j], None)

legend = axs[0][1].legend(loc='upper center', bbox_to_anchor=(0., 1.3), ncol=2, shadow=True)
for i in range(ncols):
    axs[nrows - 1][i].set_xlabel(r"$\omega [t]$")
for i in range(nrows):
    axs[i][0].set_ylabel(r"$\mathcal{A}(\omega) [t^{-1}]$")
axs[0][0].title.set_text("Square")
axs[0][1].title.set_text("Simple cubic")

axs[0][0].text(4.4, 0.55, "(a.1)\n$N_\\mathrm{\\gamma} = 300$")
axs[1][0].text(4.4, 0.55, "(b.1)\n$N_\\mathrm{\\gamma} = 1500$")
axs[2][0].text(4.4, 0.55, "(c.1)\n$N_\\mathrm{\\gamma} = 3000$")
axs[3][0].text(4.4, 0.55, "(d.1)\n$N_\\mathrm{\\gamma} = 6000$")
axs[0][1].text(7., 0.55, "(a.2)\n$N_\\mathrm{\\gamma} = 300$")
axs[1][1].text(7., 0.55, "(b.2)\n$N_\\mathrm{\\gamma} = 1500$")
axs[2][1].text(7., 0.55, "(c.2)\n$N_\\mathrm{\\gamma} = 3000$")
axs[3][1].text(7., 0.55, "(d.2)\n$N_\\mathrm{\\gamma} = 6000$")

fig.tight_layout()
plt.savefig(f"plots/{os.path.basename(__file__).split('.')[0]}.pdf")
#plt.show()