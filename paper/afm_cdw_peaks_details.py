import numpy as np
import matplotlib.pyplot as plt

import os, sys
if os.name == "nt":
    sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + r"\python")
else:
    sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + "/python")

import lib.continued_fraction as cf
from lib.iterate_containers import *
from lib.extract_key import *
import lib.resolvent_peak as rp
import lib.plot_settings as ps


Ts = np.array([0.])
Us = np.array([0.0001, 0.00015, 0.0002, 0.0003, 0.0004, 0.0005, 0.0007, 0.001, 0.0015, 0.002, 0.003, 0.004, 0.005, 0.007, 
                0.01, 0.015, 0.02, 0.03, 0.04, 0.05, 0.07, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4])
Us_cube = np.array([0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0])
Vs = np.array([1.0])


folders = ["../data/modes/square/dos_3k/", "../data/modes/cube/dos_3k/"]
element_names = ["a", "a+b", "a+ib"]

name_suffices = ["AFM", "CDW"]
nrows=2
fig, axs = plt.subplots(nrows=nrows, ncols=2, figsize=(12.8, 8.2), sharex="col", gridspec_kw=dict(hspace=0))
plotters = np.empty((nrows,2), dtype=ps.CURVEFAMILY)

for i in range(2):
    for j in range(nrows):
        plotters[j][i] = ps.CURVEFAMILY(6, axis=axs[j][i])
        plotters[j][i].set_individual_colors("nice")
        plotters[j][i].set_individual_linestyles(["-", "-"])
        plotters[j][i].set_individual_markerstyles(["v", "^"])
    if i == 0:
        u_data = 4.0 + Us
    else:
        Us = np.append(Us, Us_cube)
        u_data = 6.0 + Us
    
    weights = np.zeros(len(u_data))
    peak_positions = np.zeros(len(u_data))
    peak_positions_to_cont = np.zeros(len(u_data))
    
    for name_suffix in name_suffices:
        counter = 0
        for T, U, V in iterate_containers(Ts, u_data, Vs):
            name = f"T={T}/U={U}/V={V}"
            cont_edges = cf.continuum_edges(f"{folders[i]}{name}", f"higgs_{name_suffix}", True)
            lower = 0 if float(V) < 1 else float(V)
            upper = 13 * float(V) + 2 if 13 * float(V) + 2 < cont_edges[0] else cont_edges[0]

            peak_positions[counter], weights[counter] = rp.analyze_peak(f"{folders[i]}{name}", f"higgs_{name_suffix}", (lower, upper))
            peak_positions_to_cont[counter] = cont_edges[0] - peak_positions[counter]
            counter += 1
        
        plotters[0][i].plot(np.log(Us), np.log(peak_positions_to_cont), label=f"Data - $\\mathcal{{A}}_\\mathrm{{{name_suffix}}} (\\omega)$")
        plotters[1][i].plot(np.log(Us), weights)

axs[0][0].title.set_text("Square")
axs[0][1].title.set_text("Simple cubic")

axs[0][0].text(0.8, 0.9, "(a.1)", transform = axs[0][0].transAxes)
axs[0][1].text(0.8, 0.9, "(a.2)", transform = axs[0][1].transAxes)
axs[1][0].text(0.8, 0.9, "(b.1)", transform = axs[1][0].transAxes)
axs[1][1].text(0.8, 0.9, "(b.2)", transform = axs[1][1].transAxes)

axs[nrows-1][0].set_xlabel(r"$\ln((U - 4V) / t)$")
axs[nrows-1][1].set_xlabel(r"$\ln((U - 6V) / t)$")
axs[0][0].set_ylabel(r"$\ln((\omega_- - \omega_0) / t)$")
axs[1][0].set_ylabel(r"$w_0$")
legend = axs[0][1].legend(loc='upper center')
fig.tight_layout()

import os
plt.savefig(f"plots/{os.path.basename(__file__).split('.')[0]}.pdf")
#plt.show()
