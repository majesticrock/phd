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

# visually scale the data for visibility
SCALE = 5
params = [ [0., -2.5, -0.1], [0., -2.5, 0.1] ]

folders = ["../data/modes/square/dos_3k/", "../data/modes/cube/dos_3k/"]
nrows = 2
ncols = 2
# ax = axs[row][col]
fig, axs = plt.subplots(nrows=nrows, ncols=ncols, figsize=(12.8, 6), sharey=True, sharex="col", gridspec_kw=dict(hspace=0, wspace=0))

plotters = np.empty((nrows, ncols), dtype=ps.CURVEFAMILY)
for i in range(nrows):
    for j in range(ncols):
        axs[i][j].set_ylim(0, .99)
        plotters[i][j] = ps.CURVEFAMILY(4, axis=axs[i][j])
        plotters[i][j].set_individual_colors("nice2")
        plotters[i][j].set_individual_dashes([(1,0), (6,6), (1,0), (6,6)])

plot_lower_lim = -0.05
plot_upper_lim = 4.25

name_suffices = ["phase_SC", "higgs_SC", "CDW", "AFM"]
labels = ["Phase", "Higgs", "CDW", "AFM"]

for j, folder in enumerate(folders):
    usage_upper_lim = 2 * plot_upper_lim if j == 0 else 3 * plot_upper_lim

    for i, name in enumerate(naming_scheme_tuples(params)):
        for k, (name_suffix, label) in enumerate(zip(name_suffices, labels)):
            data, data_real, w_lin, res = cf.resolvent_data(f"{folder}{name}", name_suffix, plot_lower_lim, usage_upper_lim, 
                                                            number_of_values=5000, xp_basis=True, imaginary_offset=1e-5, messages=False)
            plotters[i][j].plot(w_lin, SCALE*data, label=label)
            
        axs[i][j].set_xlim(plot_lower_lim, usage_upper_lim)
        res.mark_continuum(axs[i][j], None)

legend = axs[0][1].legend(loc='upper center', bbox_to_anchor=(0., 1.25), ncol=2, shadow=True)
for i in range(ncols):
    axs[nrows - 1][i].set_xlabel(r"$\omega [t]$")
for i in range(nrows):
    axs[i][0].set_ylabel(r"$\mathcal{A}(\omega + i0^+) [t^{-1}]$")
axs[0][0].title.set_text("Square")
axs[0][1].title.set_text("Simple cubic")

axs[0][1].text(10.6, 0.87, r"(a) SC")
axs[1][1].text(10.6, 0.87, r"(b) CDW")

fig.tight_layout()
plt.savefig("plots/resolvent_overview_SC_CDW.pdf")