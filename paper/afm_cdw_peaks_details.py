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
from lib.ez_fit import *

Ts = np.array([0.])
Us_square = np.array([0.225, 0.25, 0.275, 0.3, 0.325, 0.35, 0.375, 0.4, 0.41, 0.42, 0.421, 0.422, 0.423]) # there is a peak at 0.424, but our numerical tools start breaking 
Us_cube   = np.array([ 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.03, 1.05, 1.06, 1.07, 1.08, 1.09, 1.1, 1.105])
# New square: 0.413, 0.415, 0.417, 0.418, 0.419
# New cube:   1.102, 1.095, 1.085
Us = np.array([Us_square, Us_cube], dtype=object)
Vs = np.array([1.0])


folders = ["../data/modes/square/dos_3k/", "../data/modes/cube/dos_3k/"]
element_names = ["a", "a+b", "a+ib"]

name_suffix = "AFM"
nrows=2
fig, axs = plt.subplots(nrows=nrows, ncols=2, figsize=(12.8, 6.4), sharex="col", gridspec_kw=dict(hspace=0))
plotters = np.empty((nrows,2), dtype=ps.CURVEFAMILY)

peak_fit_params = [ {"range": 1e-6, "begin_offset": 5e-10, "imaginary_offset": 1e-6},
                    {"range": 5e-8, "begin_offset": 1e-12, "imaginary_offset": 1e-8}]

for i in range(2):
    for j in range(nrows):
        plotters[j][i] = ps.CURVEFAMILY(6, axis=axs[j][i])
        plotters[j][i].set_individual_colors("nice")
        plotters[j][i].set_individual_linestyles(["-", ""])
        plotters[j][i].set_individual_markerstyles(["", "^"])
    if i == 0:
        u_data = 4.0 - Us[0]
    else:
        u_data = 6.0 - Us[1]
    
    weights = np.zeros(len(u_data))
    peak_positions = np.zeros(len(u_data))
    peak_positions_to_cont = np.zeros(len(u_data))
    counter = 0
    for T, U, V in iterate_containers(Ts, u_data, Vs):
        name = f"T={T}/U={round(U, 4)}/V={V}"
        cont_edges = cf.continuum_edges(f"{folders[i]}{name}", f"higgs_{name_suffix}", True)
        lower = 0 if float(V) < 1 else float(V)
        upper = 13 * float(V) + 2 if 13 * float(V) + 2 < cont_edges[0] else cont_edges[0]

        peak_positions[counter], weights[counter] = rp.analyze_peak(f"{folders[i]}{name}", f"higgs_{name_suffix}", (lower, upper), reversed=True, **(peak_fit_params[0 if counter < len(u_data) - 5 else 1]))
        peak_positions_to_cont[counter] = cont_edges[0] - peak_positions[counter]
        counter += 1
    
    u_log = np.log(u_data - (3.575 if i == 0 else 4.89))
    peak_positions_to_cont = np.log(peak_positions_to_cont)
    
    popt1, pcov1 = ez_linear_fit(u_log, peak_positions_to_cont, plotters[0][i], label="Fit")
    popt2, pcov2 = ez_linear_fit(u_log, weights, plotters[1][i], label="Fit")
    
    axs[0][i].text(0.05, 0.8, f"a = {popt1[0]:.4f}", transform = axs[0][i].transAxes)
    axs[1][i].text(0.05, 0.8, f"a = {popt2[0]:.4f}", transform = axs[1][i].transAxes)
    axs[0][i].text(0.05, 0.7, f"b = {popt1[1]:.4f}", transform = axs[0][i].transAxes)
    axs[1][i].text(0.05, 0.7, f"b = {popt2[1]:.4f}", transform = axs[1][i].transAxes)
    
    plotters[0][i].plot(u_log, peak_positions_to_cont, label=f"Data")
    plotters[1][i].plot(u_log, weights)

axs[0][0].title.set_text("Square")
axs[0][1].title.set_text("Simple cubic")

axs[0][0].text(0.89, 0.75, "(a.1)", transform = axs[0][0].transAxes)
axs[0][1].text(0.89, 0.75, "(a.2)", transform = axs[0][1].transAxes)
axs[1][0].text(0.89, 0.75, "(b.1)", transform = axs[1][0].transAxes)
axs[1][1].text(0.89, 0.75, "(b.2)", transform = axs[1][1].transAxes)

axs[nrows-1][0].set_xlabel(r"$\ln((4V - U) / t)$")
axs[nrows-1][1].set_xlabel(r"$\ln((6V - U) / t)$")
axs[0][0].set_ylabel(r"$(\omega_- - \omega_0) / t$")
axs[1][0].set_ylabel(r"$W_0$")
legend = axs[0][1].legend(loc='lower right')
fig.tight_layout()

plt.savefig(f"plots/{os.path.basename(__file__).split('.')[0]}.pdf")
#plt.show()
