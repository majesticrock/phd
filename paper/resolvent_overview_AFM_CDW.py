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

prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']

params = [ [0., 4.1, 1.], [0., 3.9, 1.] ]

use_XP = True

folders = ["../data/modes/square/dos_3k/", "../data/modes/cube/dos_3k/"]
nrows = 2
ncols = 2
# ax = axs[row][col]
fig, axs = plt.subplots(nrows=nrows, ncols=ncols, figsize=(12.8, 6), sharey=True, sharex="col", gridspec_kw=dict(hspace=0, wspace=0))

plotters = np.empty((nrows, ncols), dtype=ps.CURVEFAMILY)
for i in range(nrows):
    for j in range(ncols):
        axs[i][j].set_ylim(0, 0.75)
        plotters[i][j] = ps.CURVEFAMILY(4, axis=axs[i][j])
        plotters[i][j].set_individual_colors("nice")
        plotters[i][j].set_individual_linestyles(["-", "--", "-", "-"])

plot_lower_lim = 0#-0.05
plot_upper_lim = 5.2

name_suffices = ["phase_SC", "higgs_SC", "CDW", "AFM"]
labels = ["Phase", "Higgs", "CDW", "AFM"]

for j, folder in enumerate(folders):
    usage_upper_lim = 2 * plot_upper_lim if j == 0 else 2.7 * plot_upper_lim
    if j == 1: 
        # the sc lattice uses a different V for the CDW-AFM border due to the changed coordination number
        params[1][2] = 0.5
        params[0][2] = 0.5
        
    for i, name in enumerate(naming_scheme_tuples(params)):
        for name_suffix, label in zip(name_suffices, labels):
            data, data_real, w_lin, res = cf.resolvent_data(f"{folder}{name}", name_suffix, plot_lower_lim, usage_upper_lim, 
                                                            number_of_values=20000, xp_basis=use_XP, imaginary_offset=1e-6, messages=False)
            plotters[i][j].plot(w_lin, data, label=label)
            
            cont = np.sqrt(res.roots[0])
            create_zoom(axs[i][j], 0.5, 0.05, 0.2, 0.9, (cont - 0.025, cont), ylim=(0, 0.35))
            
        axs[i][j].set_xlim(plot_lower_lim, usage_upper_lim)
        res.mark_continuum(axs[i][j])


legend = axs[0][0].legend(loc="upper center")
for i in range(ncols):
    axs[nrows - 1][i].set_xlabel(r"$z / t$")
for i in range(nrows):
    axs[i][0].set_ylabel(r"$\mathcal{A}(\omega)$ / a.u.")
axs[0][0].title.set_text("Square lattice")
axs[0][1].title.set_text("Simple cubic lattice")

axs[0][1].text(12.5, 0.65, r"(a)")
axs[1][1].text(12.5, 0.65, r"(b)")

fig.tight_layout()
plt.savefig("plots/resolvent_overview_AFM_CDW.pdf")