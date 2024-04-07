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

poster_plot = False
if poster_plot:
    legend_args = { "bbox_to_anchor" : (-0.05, 1.19), "columnspacing" : 1, "fontsize" : 22 }
else:
    legend_args = { "bbox_to_anchor" : (0., 1.15)}

params = [ [0., -2.5, -0.1], [0., -2.5, 0.1], [0., -2.5, 0.5] ]
folders = ["../data/modes/square/dos_6000/", "../data/modes/cube/dos_6000/"]
nrows = 3
ncols = 2
# ax = axs[row][col]
fig, axs = plt.subplots(nrows=nrows, ncols=ncols, figsize=(12.8, 8), sharey=True, sharex="col", gridspec_kw=dict(hspace=0, wspace=0))

plotters = np.empty((nrows, ncols), dtype=ps.CURVEFAMILY)
for i in range(nrows):
    for j in range(ncols):
        axs[i][j].set_ylim(0, .99)
        plotters[i][j] = ps.CURVEFAMILY(5, axis=axs[i][j])
        plotters[i][j].set_individual_colors("GPT")
        plotters[i][j].set_individual_dashes([ [1,0], [1.5,1.5], [3, 1.5], [5,2,1,2], [7,3,2,2,3,2] ])

name_suffices = ["phase_SC", "higgs_SC", "CDW", "AFM", "higgs_AFM_trans"]
labels = ["Phase", "Higgs", "CDW", "l.AFM", "t.AFM"]

plot_lower_lim = -0.05
plot_upper_lim = 4.25

for j, folder in enumerate(folders):
    usage_upper_lim = 2 * plot_upper_lim if j == 0 else 3 * plot_upper_lim

    for i, name in enumerate(naming_scheme_tuples(params)):
        resolvents = np.empty(len(name_suffices), dtype=cf.ContinuedFraction)
        for k, (name_suffix, label) in enumerate(zip(name_suffices, labels)):
            data, data_real, w_lin, resolvents[k] = cf.resolvent_data(f"{folder}{name}", name_suffix, plot_lower_lim, usage_upper_lim, 
                                                            number_of_values=7000, xp_basis=True, imaginary_offset=1e-5, messages=False)
            plotters[i][j].plot_with_peak(w_lin, data, label=r"$\mathcal{A}_\mathrm{" + label + r"} (\omega)$")
            
        axs[i][j].set_xlim(plot_lower_lim, usage_upper_lim)
        if i == 2:
            cont = np.sqrt(resolvents[0].roots[0])
            zoomed_xlim = (cont - 0.08, cont + 0.03) 
            axins = create_zoom(axs[i][j], 0.1, 0.29, 0.275, 0.66, zoomed_xlim, ylim=(0, 0.55), 
                                y_funcs=[lambda x, res=res: res.spectral_density(x + 1e-5j) for res in resolvents],
                                skip_lines=[1, 3, 5, 7, 9], yticks=[0, 0.2, 0.4], mark_inset=False)
            resolvents[0].mark_continuum(axins, None)
            
        resolvents[0].mark_continuum(axs[i][j], None)
        

legend = axs[0][1].legend(loc='upper center', ncol=5, shadow=True, **legend_args)

for i in range(ncols):
    axs[nrows - 1][i].set_xlabel(r"$\omega [t]$")
for i in range(nrows):
    axs[i][0].set_ylabel(r"$\mathcal{A}(\omega) [t^{-1}]$")
axs[0][0].set_title("Square lattice", pad=22 if not poster_plot else 25)
axs[0][1].set_title("Simple cubic lattice", pad=22 if not poster_plot else 25)

axs[0][0].text(0.93, 0.59, "(a.1) SC\n$V=-0.1t$", transform = axs[0][0].transAxes, ma="right", ha="right")
axs[1][0].text(0.93, 0.59, "(b.1) CDW\n$V=0.1t$", transform = axs[1][0].transAxes, ma="right", ha="right")
axs[2][0].text(0.93, 0.59, "(c.1) CDW\n$V=0.5t$", transform = axs[2][0].transAxes, ma="right", ha="right")
axs[0][1].text(0.93, 0.59, "(a.2) SC\n$V=-0.1t$", transform = axs[0][1].transAxes, ma="right", ha="right")
axs[1][1].text(0.93, 0.59, "(b.2) CDW\n$V=0.1t$", transform = axs[1][1].transAxes, ma="right", ha="right")
axs[2][1].text(0.93, 0.59, "(c.2) CDW\n$V=0.5t$", transform = axs[2][1].transAxes, ma="right", ha="right")

fig.tight_layout()
plt.savefig(f"plots/{os.path.basename(__file__).split('.')[0]}.pdf")
plt.savefig(f"plots/{os.path.basename(__file__).split('.')[0]}.svg")
#plt.show()
