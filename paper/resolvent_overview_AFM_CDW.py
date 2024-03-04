from pydoc import resolve
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

params = [  [ [0., 4.85, 1.2], [0., 4.75, 1.2], [0., 6.2, 1.2], [0., 3.4, 1.2] ],
            [ [0., 4.85, 0.8], [0., 4.75, 0.8], [0., 6.2, 0.8], [0., 3.4, 0.8] ]]

folders = ["../data/modes/square/dos_6000/", "../data/modes/cube/dos_6000/"]
nrows = 4
ncols = 2
# ax = axs[row][col]
fig, axs = plt.subplots(nrows=nrows, ncols=ncols, figsize=(12.8, 10), sharey=True, sharex="col", gridspec_kw=dict(hspace=0, wspace=0))

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
plot_upper_lim = 5.35

for j, folder in enumerate(folders):
    usage_upper_lim = 2 * plot_upper_lim if j == 0 else 2.8 * plot_upper_lim
        
    for i, name in enumerate(naming_scheme_tuples(params[j])):
        if i >= nrows: continue
        resolvents = np.empty(len(name_suffices), dtype=cf.ContinuedFraction)
        for k, (name_suffix, label) in enumerate(zip(name_suffices, labels)):
            data, data_real, w_lin, resolvents[k] = cf.resolvent_data(f"{folder}{name}", name_suffix, plot_lower_lim, usage_upper_lim, 
                                                            number_of_values=10000, xp_basis=True, imaginary_offset=1e-5, messages=False)
            plotters[i][j].plot_with_peak(w_lin, data, label=r"$\mathcal{A}_\mathrm{" + label + r"} (\omega)$")
            
        cont = np.sqrt(resolvents[0].roots[0])
        if i < 2:
            zoomed_xlim = np.array([cont - 0.035, cont + 0.025]) if j == 1 else np.array([cont - 0.09, cont + 0.05]) 
        else:
            zoomed_xlim = np.array([cont - 0.07, cont + 0.025]) if j == 1 else np.array([cont - 0.14, cont + 0.05]) 
            
        zoom_xpos = 0.4 if i < 2 or j == 1 else 0.1
        axins = create_zoom(axs[i][j], zoom_xpos, 0.29, 0.275, 0.66, zoomed_xlim, ylim=(0, 0.55), 
                            y_funcs=[lambda x, res=res: res.spectral_density(x + 1e-5j) for res in resolvents],
                            skip_lines=[1, 3, 5, 7, 9], yticks=[0, 0.2, 0.4], mark_inset=False)
        axs[i][j].set_xlim(plot_lower_lim, usage_upper_lim)
        resolvents[0].mark_continuum(axs[i][j], None)
        resolvents[0].mark_continuum(axins, None)


legend = axs[0][1].legend(loc='upper center', bbox_to_anchor=(0., 1.3), ncol=2, shadow=True)

legend = axs[0][1].legend(loc='upper center', bbox_to_anchor=(0., 1.15), ncol=5, shadow=True)
for i in range(ncols):
    axs[nrows - 1][i].set_xlabel(r"$\omega [t]$")
for i in range(nrows):
    axs[i][0].set_ylabel(r"$\mathcal{A}(\omega) [t^{-1}]$")
axs[0][0].set_title("Square lattice", pad=22)
axs[0][1].set_title("Simple cubic lattice", pad=22)

axs[0][0].text(7.6, 0.6, "(a.1) AFM\n$U = 4.85t$")
axs[0][1].text(10.6, 0.6, "(a.2) AFM\n$U = 6.2t$")

axs[1][1].text(10.6, 0.6, "(b.2) CDW\n$U = 3.4t$")
axs[1][0].text(7.6, 0.6, "(b.1) CDW\n$U = 4.75t$")

axs[2][0].text(7.6, 0.6, "(c.1) AFM\n$U = 4.85t$")
axs[2][1].text(10.6, 0.6, "(c.2) AFM\n$U = 6.2t$")

if nrows > 3: axs[3][0].text(7.6, 0.6, "(d.1) CDW\n$U = 4.75t$")
if nrows > 3: axs[3][1].text(10.6, 0.6, "(d.2) CDW\n$U = 3.4t$")

fig.tight_layout()
plt.savefig(f"plots/{os.path.basename(__file__).split('.')[0]}.pdf")
#plt.show()
